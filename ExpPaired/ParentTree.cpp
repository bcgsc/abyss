#include "ParentTree.h"
#include "PackedSeq.h"
#include "Timer.h"

ParentTree::ParentTree(const PackedSeq& rootNode, extDirection dir, ISequenceCollection* pSC, int maxDepth) : m_root(rootNode), m_pSC(pSC), m_dir(dir), m_maxDepth(maxDepth)
{
	Timer timer("tree build");
	// Construct the tree by exploring out from the root node and adding children up to max depth
	explore(rootNode, 0);
	
	// Now that the tree is constructed, cache the relationship from a node to a root-child (second level node)
	// most algorithms will operate on this table
	cacheNodeToRootChild();
}

void ParentTree::addScore(const PackedSeq& node, double weight)
{
	// Check if the seq or its RC should be added
	PackedSeq startNode;
	bool found = false;
	RelationshipTable::iterator parents = m_childToParent.find(node);
	if(parents != m_childToParent.end())
	{
		startNode = node;
		found = true;
	}
	else
	{
		PackedSeq rc = reverseComplement(node);
		RelationshipTable::iterator parentsRC = m_childToParent.find(rc);
		if(parentsRC != m_childToParent.end())
		{
			startNode = rc;
			found = true;
		}	
	}	
	
	if(found)
	{
		// generate all the sequences to update
		PSeqSet updateSet;
		generateReachableSet(startNode, updateSet);
		
		// update all the sequences
		for(PSeqSet::iterator iter = updateSet.begin(); iter != updateSet.end(); ++iter)
		{
			bool actualseq = false;
			
			if(*iter == startNode)
			{
				actualseq = true;
			}
			
			m_nodeScores.addScore(*iter, actualseq, weight);
		}
	}
}

TreeScore ParentTree::getScore(const PackedSeq& node)
{
	TreeScore score = m_nodeScores.getScore(node);
	return score;
}

void ParentTree::explore(const PackedSeq& parent, int depth)
{
	
	m_nodeScores.createNode(parent);
		
	if(depth >= m_maxDepth)
	{
		return;
	}

	// Find the children of this sequence
	ExtensionRecord extRecord;
	int multiplicity;
	bool seqFound = m_pSC->getSeqData(parent, extRecord, multiplicity);
	if(seqFound)
	{	

		PSequenceVector children;
		AssemblyAlgorithms::generateSequencesFromExtension(parent, m_dir, extRecord.dir[m_dir], children);

		for(PSequenceVector::iterator iter = children.begin(); iter != children.end(); iter++)
		{
			bool shouldRecurse = (m_childToParent.find(*iter) == m_childToParent.end()) ? true : false;
			
			// Add a child->parent record to the table
			addToTable(m_childToParent, *iter, parent);
			addToTable(m_parentToChild, parent, *iter);

			// recurse if the iter is not already in the table
			if(shouldRecurse)
			{
				explore(*iter, depth + 1);
			}
		}
	}
	else
	{
		printf("seq not found! dead end\n");
	}
}

//
// Cache the relationship between a node and a root-child
//
void ParentTree::cacheNodeToRootChild()
{
	Timer t("cacheNodeToRC");
	PSeqSet rootChildren = getRootChildren();
	
	PSeqQueue nodeQueue;
		
	// base case, make the root node children point to themselves
	for(PSeqSet::iterator rcIter = rootChildren.begin(); rcIter != rootChildren.end(); ++rcIter)
	{
		addToTable(m_nodeToRootChild, *rcIter, *rcIter);
		addChildrenToQueue(*rcIter, nodeQueue);
	}
	
	PSeqSet loopGuard;
	
	// Perform a breadth-first search, setting each node to have the values of the sum of its parents
	while(!nodeQueue.empty())
	{
		const PackedSeq currSeq = nodeQueue.front();
		
		// get the list of this node's parents
		RelationshipTable::iterator parents = m_childToParent.find(currSeq);
		if(parents != m_childToParent.end())
		{
			for(PSeqSet::iterator pIter = parents->second.begin(); pIter != parents->second.end(); ++pIter)
			{
				RelationshipTable::iterator parentNodeData = m_nodeToRootChild.find(*pIter);
				if(parentNodeData != m_nodeToRootChild.end())
				{
					for(PSeqSet::iterator pnIter = parentNodeData->second.begin(); pnIter != parentNodeData->second.end(); ++pnIter)
					{
						addToTable(m_nodeToRootChild, currSeq, *pnIter);
					}
				}
			}
		}
		else
		{
			assert(false);	
		}
		
		// Add this node's children to the queue if it has not been seen before
		if(loopGuard.find(currSeq) == loopGuard.end())
		{
			addChildrenToQueue(currSeq, nodeQueue);
			loopGuard.insert(currSeq);	
		}
		
		nodeQueue.pop();
	}
}

// 
// Score the root's children only by using the cache table
//
void ParentTree::scoreRootChildrenOnly(const PackedSeq& node, double weight)
{
	scoreRootChildrenOnlyInternal(node, weight);
	scoreRootChildrenOnlyInternal(reverseComplement(node), weight);
}


void ParentTree::scoreRootChildrenOnlyInternal(const PackedSeq& node, double weight)
{
	RelationshipTable::iterator rc = m_nodeToRootChild.find(node);
	if(rc != m_nodeToRootChild.end())
	{
		for(PSeqSet::iterator iter = rc->second.begin(); iter != rc->second.end(); iter++)
		{
			bool actualseq = false;
			if(*iter == node)
			{
				actualseq = true;
			}
			
			m_nodeScores.addScore(*iter, actualseq, weight);
		}
	}	
}

//
// Clear the scores of the root's children
//
void ParentTree::clearRootChildrenScores()
{
	RelationshipTable::iterator children = m_parentToChild.find(m_root);
	if(children != m_parentToChild.end())
	{
		for(PSeqSet::iterator iter = children->second.begin(); iter != children->second.end(); iter++)
		{
			m_nodeScores.clearScore(*iter);
		}
	}
}

//
// build a path backwards through the tree
//
void ParentTree::buildBackwards(const PackedSeq& begin, const PackedSeq& end, PSequenceVector& currSeqs, std::vector<PSequenceVector>& paths, size_t limit, size_t maxPaths)
{
	if(begin == end)
	{
		// This path is good, add it
		paths.push_back(currSeqs);
		return;
	}
	
	if(currSeqs.size() > limit)
	{
		return;	
	}
	
	if(paths.size() > maxPaths)
	{
		return;
	}
	
	// get the parents of begin
	RelationshipTable::iterator parents = m_childToParent.find(begin);
	if(parents != m_childToParent.end())
	{	
		if(parents->second.size() == 1)
		{
			// tail recursive case
			currSeqs.push_back(*parents->second.begin());
			return buildBackwards(*parents->second.begin(), end, currSeqs, paths, limit, maxPaths);
		}
		else
		{
			// fork
			for(PSeqSet::iterator iter = parents->second.begin(); iter != parents->second.end(); iter++)
			{
				PSequenceVector forkSeqs = currSeqs;
				forkSeqs.push_back(*iter);
				buildBackwards(*iter, end, forkSeqs, paths, limit, maxPaths);
			}
		}
	}
	else
	{
		return;
	}
}

//
//
//
void ParentTree::generateReachableSet(const PackedSeq& node, PSeqSet& updateSet)
{
	// add the current node to the set
	updateSet.insert(node);
	
	// add the parents of the node
	RelationshipTable::iterator parents = m_childToParent.find(node);
	if(parents != m_childToParent.end())
	{
		for(PSeqSet::iterator iter = parents->second.begin(); iter != parents->second.end(); iter++)
		{
			if(updateSet.find(*iter) == updateSet.end())
			{
				generateReachableSet(*iter, updateSet);
			}
		}
	}	
}

void ParentTree::resetScore()
{
	m_nodeScores.reset();
}

void ParentTree::addToTable(RelationshipTable& table, const PackedSeq& key, const PackedSeq& value)
{
	table[key].insert(value);
}

PSeqSet ParentTree::getRootChildren()
{
	RelationshipTable::iterator children = m_parentToChild.find(m_root);
	if(children != m_parentToChild.end())
	{
		return children->second;
	}
	else
	{
		PSeqSet dummy;
		return dummy;	
	}
}

//
//
//
void ParentTree::addChildrenToQueue(const PackedSeq& node, PSeqQueue& queue)
{
	// add its children to the queue
	RelationshipTable::iterator children = m_parentToChild.find(node);
	if(children != m_parentToChild.end())
	{
		for(PSeqSet::iterator cIter = children->second.begin(); cIter != children->second.end(); cIter++)
		{
			queue.push(*cIter);
		}
	}	
}

void ParentTree::printRootChildren()
{
	printf("The root node is: %s\n", m_root.decode().c_str());
	printf("Its children:\n");
	
	RelationshipTable::iterator children = m_parentToChild.find(m_root);
	if(children != m_parentToChild.end())
	{
		for(PSeqSet::iterator iter = children->second.begin(); iter != children->second.end(); iter++)
		{
			TreeScore score = m_nodeScores.getScore(*iter);
			printf("\t%s score: (%lf,%lf)\n", iter->decode().c_str(), score.indScore, score.subtreeScore);
		}
	}		
}

void ParentTree::print(int maxDepth)
{
	printf("The tree has %zu (%zu) nodes\n", m_childToParent.size(), m_parentToChild.size());
	printf("The root node is: %s\n", m_root.decode().c_str());
	printf("Root has %zu children\n", m_parentToChild.find(m_root)->second.size());
	printRecursive(m_root, 1, maxDepth);	
}

void ParentTree::printRecursive(const PackedSeq& node, int depth, int maxDepth)
{
	if(depth > maxDepth)
	{
		return;
	}
	
	for(int i = 0; i < depth; i++)
	{
		printf(" ");
	}
	
	TreeScore score = m_nodeScores.getScore(node);
	
	
	RelationshipTable::iterator children = m_parentToChild.find(node);
	if(children != m_parentToChild.end())
	{
		printf("%s (%lf, %lf)-> (%zu)\n", node.decode().c_str(), score.indScore, score.subtreeScore, children->second.size());
		for(PSeqSet::iterator iter = children->second.begin(); iter != children->second.end(); iter++)
		{
			printRecursive(*iter, depth+1, maxDepth);
		}
	}
	else
	{
		printf("%s (%lf, %lf)\n", node.decode().c_str(), score.indScore, score.subtreeScore);
	}
}



