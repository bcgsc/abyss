#include "ParentTree.h"
#include "PackedSeq.h"
#include "Timer.h"

ParentTree::ParentTree(const PackedSeq& rootNode, extDirection dir, ISequenceCollection* pSC, int maxDepth) : m_root(rootNode), m_pSC(pSC), m_dir(dir), m_maxDepth(maxDepth)
{
	Timer timer("tree build");
	// Construct the tree by exploring out from the root node and adding children up to max depth
	explore(rootNode, 0);
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
		generateUpdateSet(startNode, updateSet);
		
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
	bool seqFound = m_pSC->getExtensions(parent, extRecord);
	if(seqFound)
	{	

		PSequenceVector children;
		AssemblyAlgorithms::generateSequencesFromExtension(parent, m_dir, extRecord.dir[m_dir], children);

		for(PSequenceVector::iterator iter = children.begin(); iter != children.end(); iter++)
		{
			// Add a child->parent record to the table
			addToTable(m_childToParent, *iter, parent);
			addToTable(m_parentToChild, parent, *iter);

				
			// recurse
			explore(*iter, depth + 1);
		}
	}
	else
	{
		printf("seq not found! dead end\n");
	}
}

//
//
//
void ParentTree::generateUpdateSet(const PackedSeq& node, PSeqSet& updateSet)
{
	// add the current node to the set
	updateSet.insert(node);
	
	// add the parents of the node
	RelationshipTable::iterator parents = m_childToParent.find(node);
	if(parents != m_childToParent.end())
	{
		for(PSeqSet::iterator iter = parents->second.begin(); iter != parents->second.end(); iter++)
		{
			generateUpdateSet(*iter, updateSet);
		}
	}	
}

void ParentTree::addToTable(RelationshipTable& table, const PackedSeq& key, const PackedSeq& value)
{
	table[key].insert(value);
}

PSeqSet ParentTree::getRootChildren()
{
	RelationshipTable::iterator children = m_parentToChild.find(m_root);
	assert(children != m_parentToChild.end());
	return children->second;
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



