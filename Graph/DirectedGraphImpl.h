#include "DirectedGraph.h"
#include "ContigNode.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib> // for exit
#include <iostream>

namespace opt {
	extern unsigned k;
};

struct SimpleDataCost
{
	/** Return the length of the specified node in k-mer. */
	size_t cost(const SimpleContigData& data)
	{
		return data.length - opt::k + 1;
	}
};
static SimpleDataCost costFunctor;

template<typename K, typename D>
void Vertex<K,D>::addEdge(VertexType* pNode, extDirection dir, bool reverse) 
{
	EdgeData data;
	data.pVertex = pNode;
	data.reverse = reverse;
	
	// Check if this edge already exists
	for(typename EdgeCollection::const_iterator edgeIter = m_edges[dir].begin(); edgeIter != m_edges[dir].end(); ++edgeIter)
	{
		if(*edgeIter == data)
		{
			std::cout << "Warning, attempted to add duplicate edge (" << m_key << ") -> " << pNode->m_key << "," << dir << "," << reverse << "\n";
			return;
		}
	}
	m_edges[dir].push_back(data);
}

//
//
//
template<typename K, typename D>
void Vertex<K,D>::removeEdge(VertexType* pNode, extDirection dir, bool reverse)
{
	// Remove the node
	bool found;
	typename EdgeCollection::iterator edgeIter = getEdge(pNode, dir, reverse, found);

	// Make sure the edge was found
	assert(edgeIter != m_edges[dir].end() && found);
	 m_edges[dir].erase(edgeIter); // slow, resizes vector but memory efficiency is most important
}
 
//
//
//
template<typename K, typename D>
typename Vertex<K,D>::EdgeCollectionIter Vertex<K,D>::getEdge(VertexType* pNode, extDirection dir, bool reverse, bool& found)
{
	EdgeData data;
	data.pVertex = pNode;
	data.reverse = reverse;
		
	// Search for the edge in the vertex collection
	EdgeCollection& currEdgeSet = m_edges[dir];
	typename EdgeCollection::iterator edgeIter = currEdgeSet.begin();
	while(edgeIter != currEdgeSet.end())
	{
		if(*edgeIter == data)
		{
			// Edge found;
			break;
		}
		++edgeIter;
	}

	found = (edgeIter != currEdgeSet.end());
	return edgeIter;
}

//
// Returns true if the described edge is the only edge in the direction
//
template<typename K, typename D>
bool Vertex<K,D>::isEdgeUnique(VertexType* pNode, extDirection dir, bool reverse)
{
	// first, make sure the edge is actually in the collection
	bool found;
	typename Vertex<K,D>::EdgeCollectionIter iter = getEdge(pNode, dir, reverse, found);
	assert(found);
	
	// if the edge is found and there is only one edge in the direction, it has to be unique
	return (numEdges(dir) == 1);
}

//
//
//
template<typename K, typename D>
EdgeDescription Vertex<K,D>::findUniqueEdge(const K& key)
{
	bool found = false;
	EdgeDescription ret;
	for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
	{
		const EdgeCollection& currEdgeSet = m_edges[dirIdx];
		for(EdgeCollectionConstIter edgeIter = currEdgeSet.begin(); edgeIter != currEdgeSet.end(); ++edgeIter)
		{
			if(edgeIter->pVertex->m_key == key)
			{
				assert(!found);
				found = true;
				ret.dir = (extDirection)dirIdx;
				ret.reverse = edgeIter->reverse;
			}
		}
	}
	assert(found);
	return ret;
}

template<typename K, typename D>
bool Vertex<K,D>::edgeExists(const K& key, extDirection dir, bool reverse)
{
	const EdgeCollection& currEdgeSet = m_edges[dir];
	for(EdgeCollectionConstIter edgeIter = currEdgeSet.begin(); edgeIter != currEdgeSet.end(); ++edgeIter)
	{
		if(edgeIter->pVertex->m_key == key && edgeIter->reverse == reverse)
		{
			return true;
		}
	}
	return false;
}

//
//
//
template<typename K, typename D>
bool Vertex<K,D>::detectSimpleCycle()
{
	// Check if any sense edge is identical to any antisense edge
	for(EdgeCollectionConstIter iter = m_edges[SENSE].begin(); iter != m_edges[SENSE].end(); ++iter)
	{
		if(m_edges[ANTISENSE].find(iter->first) != m_edges[ANTISENSE].end())
		{
			return true;
		}
	}
	return false;
}
	
template<typename K, typename D>
void Vertex<K,D>::printLinks(const EdgeCollection& collection) const
{
	for(EdgeCollectionConstIter iter = collection.begin(); iter != collection.end(); ++iter)
	{
		std::cout << iter->pVertex->m_key << "(" << iter->reverse << ") ";
	}
}

template<typename K, typename D>
void Vertex<K,D>::printEdges() const
{
	for(size_t idx = 0; idx <= 1; ++idx)
	{
		printf("Dir(%zu): ",  idx);
		printLinks(m_edges[idx]);
		printf("\n");
	}
}

//
// Graph Implementation
//
template<typename D>
DirectedGraph<D>::~DirectedGraph()
{
	// delete all the vertices in the table
	for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		delete *iter;
		*iter = NULL;
	}
}

template<typename D>
void DirectedGraph<D>::addEdge(
		const LinearNumKey& parent, extDirection dir,
		const ContigNode& child)
{
	VertexType* pParentVertex = findVertex(parent);
	assert(pParentVertex != NULL);
	VertexType* pChildVertex = findVertex(child.id());
	assert(pChildVertex != NULL);
	pParentVertex->addEdge(pChildVertex, dir, child.sense());
}

template<typename D>
typename DirectedGraph<D>::VertexType* DirectedGraph<D>::addVertex(const LinearNumKey& key, const D& data)
{
	//assert(findVertex(key) == NULL);
	
	// Create new vertex
	VertexType* ptr = new VertexType(key, data);
	
	// Insert into the table
	
	// ensure contiguitiy of the keys
	assert(m_vertexTable.size() == key);
	m_vertexTable.push_back(ptr);
	
	return ptr;
}

template<typename D>
typename DirectedGraph<D>::VertexType* DirectedGraph<D>::findVertex(const LinearNumKey& key) const
{
	// Look up the key in the vertex table
	if(key >= m_vertexTable.size())
	{
		std::cerr << "Vertex key is invalid in find! (" << key << ")\n";
		assert(false);
	}
	
	return m_vertexTable[key];
}

// Count all the edges in all the nodes
template<typename D>
size_t DirectedGraph<D>::countEdges() const
{
	size_t sum = 0;
	for(VertexTableConstIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		sum += (*iter)->countEdges();
	}
	return sum;
}

// Print out a node
template<typename D>
void DirectedGraph<D>::printVertex(const LinearNumKey& key) const
{
	// find the vertex
	VertexType* pVertex = findVertex(key);
	
	// should not be null
	assert(pVertex != NULL);
	
	// Print the vertex's key
	std::cout << "Vertex key: " << pVertex->m_key << std::endl;
	pVertex->printEdges();
}

//
// Attempt to reduce the data set using paired reads
//
template<typename D>
template<class ResolveFunctor>
size_t DirectedGraph<D>::reducePaired(ResolveFunctor& resolver)
{
	size_t numMerged = 0;
	
	for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		//if(iter->second->m_data.m_copyNumber == 1)
		if(iter->second->m_data.getLength() > 500)
		{
			bool merged = false;
			bool stop = false;
			while(!stop)
			{
				if(resolver.resolve(this, iter->first))
				{
					merged = true;
				}
				else
				{
					stop = true;
				}
			}
			
			if(!merged)
			{
				printf("***** COULD NOT MERGE %s *****\n", iter->first.c_str());
			}
		}
	}

	return numMerged;
		
}

//
// Remove transitivity in the data set by merging and 
//
template<typename D>
template<class Functor>
size_t DirectedGraph<D>::removeTransitivity(Functor dataMerger)
{

	size_t numMerged = 0;

	for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		VertexType* pVertex = iter->second;
		// Check if this node has any transitive edges
		for(size_t idx = 0; idx < NUM_DIRECTIONS; ++idx)
		{
			typename VertexType::EdgeCollection currEdges = pVertex->m_edges[idx];
			
			// Check if this direction can be merged
			if(currEdges.size() == 1)
			{
				//printf("attempting merge for %s\n", pVertex->m_key.c_str());
				//pVertex->printEdges();
				
				// Get the vertex to merge with
				
				// This statement is only valid because size == 1
				VertexType* pPartner = currEdges.begin()->pVertex;
				bool parentRev = currEdges.begin()->reverse;
				
				extDirection parentDir = (extDirection)idx;
				
				// Get the direction from the child back to the parent
				extDirection childDir = EdgeDescription::getTwinDir(parentDir, parentRev);
				
				// remove the child if the edge back to the parent is unique
				// This implies that the parent has a single extension to the child and the child
				// has a single extension to the parent so after the append the child will be redundant
				bool removeChild = pPartner->isEdgeUnique(pVertex, childDir, parentRev);
				
				// sanity check
				if(pVertex == pPartner)
				{
					continue;
				}
				
				// attempt the merge
				if(merge(pVertex, pPartner, (extDirection)idx, parentRev, removeChild, dataMerger))
				{
					numMerged++;
				}
			}
		}		
	}
	
	return numMerged;
}

//
//
//
template<typename D>
template<class Functor>
bool DirectedGraph<D>::mergeWrapper(const LinearNumKey& key1, const LinearNumKey& key2, bool forceRemove, Functor dataMerger)
{
	VertexType* pParent = findVertex(key1);
	VertexType* pChild = findVertex(key2);

	EdgeDescription edgeDesc = pParent->findUniqueEdge(key2);
	
	extDirection childDir = EdgeDescription::getTwinDir(edgeDesc.dir, edgeDesc.reverse);
	
	// Should the child be removed?
	bool removeChild = false;
	
	// is there only one edge between the parent and child?
	if(pChild->numEdges(childDir) == 1 || forceRemove)
	{
		removeChild = true;
	}
	
	return merge(pParent, pChild, edgeDesc.dir, edgeDesc.reverse, removeChild, dataMerger);
}

//
// Append a copy of the child vertex into the parent and update all the links accordingly
//
template<typename D>
template<class Functor>
bool DirectedGraph<D>::merge(VertexType* pParent, VertexType* pChild, const extDirection parentsDir, const bool parentsReverse, bool removeChild, bool usableChild, Functor dataMerger)
{
	// Determine the relative orientation of the nodes
	// If they point towards each other, the child vertex must be flipped around
	// ---SENSE---> <---SENSE = FLIP
	// ----SENSE---> <---ANTISENSE = OK
	
	/*
	printf("trying to append %s to %s (%d %d)\n", pParent->m_key.c_str(), pChild->m_key.c_str(), parentsDir, parentsReverse);
	printf("pre merge\n");
	printVertex(pParent->m_key);
	printVertex(pChild->m_key);
	*/
	
	LinearNumKey parentKey = pParent->m_key;
	LinearNumKey childKey = pChild->m_key;
	
	/*
	// check if this node has been merged before
	if(pParent->m_mergeRecord.find(childKey) != pParent->m_mergeRecord.end())
	{
		// TODO: loop found, handle this better
		assert(false);
		//do not merge, loop
		return false;
	}
	*/
	// Get the actual edge of the parent
	//bool edgeFound;
	//typename VertexType::EdgeCollectionIter parentEdge = pParent->getEdge(pChild, parentsDir, parentsReverse, edgeFound);
	//assert(edgeFound);
	
	// Compute the direction the child's edge SHOULD be in
	extDirection expectedChildsDir = (parentsReverse) ? parentsDir : !parentsDir;
	bool expectedChildsReverse = parentsReverse;
	
	// Get the child's edge
	bool edgeFound;
	typename VertexType::EdgeCollectionIter childEdge = pChild->getEdge(pParent, expectedChildsDir, expectedChildsReverse, edgeFound);
	assert(edgeFound);

	// merge the data using the functor object
	dataMerger.merge(pParent->m_key, pParent->m_data, pChild->m_key, pChild->m_data, (extDirection)parentsDir, parentsReverse, removeChild, usableChild);
	
	// update all the edges affected
	
	// set the parents edges to the edges to that of the child, taking care to set the correct reverse flag and directionality
	
	// as this link is now considered to be resolved, remove the link to the parent from all its children in this direction
	typename VertexType::EdgeCollection& parentsEdges = pParent->m_edges[parentsDir];
	for(typename VertexType::EdgeCollectionIter peIter = parentsEdges.begin(); peIter != parentsEdges.end(); ++peIter)
	{
		// set the expected direction
		extDirection expectedDir = (peIter->reverse) ? parentsDir : !parentsDir;
		bool expectedReverse = peIter->reverse;
		
		// remove the edge
		peIter->pVertex->removeEdge(pParent, expectedDir, expectedReverse);
	}
	
	// Clear the parent's edges in this direction
	pParent->m_edges[parentsDir].clear();
		
	// for each edge of the child in the opposite direction of the parent, add to the parent
	extDirection childUpdateEdgeDir = !expectedChildsDir;
	typename VertexType::EdgeCollection& childOppEdges = pChild->m_edges[childUpdateEdgeDir];
	
	for(typename VertexType::EdgeCollectionIter ceIter = childOppEdges.begin(); ceIter != childOppEdges.end(); ++ceIter)
	{
		// If the child is opposite complement of the parent, flip the reverse flag for the add
		bool newEdgeReversed =  parentsReverse != ceIter->reverse;
		//printf("Parent rev: %d Child rev: %d\n", parentEdge->second.reverse, ceIter->second.reverse);
		//printf("adding edge from %s to %s in dir %d with comp %d\n", pParent->m_key.c_str(), ceIter->first->m_key.c_str(), parentsDir, newEdgeReversed);
		pParent->addEdge(ceIter->pVertex, parentsDir, newEdgeReversed);
		
		// compute the directionality of the return edge
		
		// the reverseness is of course the same as the parent node
		// the direction of the return node is xor'd the reverseness of the parent direction
		extDirection returnEdgeDir = (newEdgeReversed) ? parentsDir : !parentsDir;

		// add the direction to the opposite vertex
		ceIter->pVertex->addEdge(pParent, returnEdgeDir, newEdgeReversed);
	}
	
	// check if the child should be removed
	if(removeChild)
	{
		// remove the vertex and update all the links to point to the merged node
		removeVertex(pChild);
	}
	
	/*
	printf("post merge\n");
	printVertex(pParent->m_key);
	//printVertex(pPartner->m_key);	
	*/
	
	// merged, return true
	return true;
}

template<typename D>
void DirectedGraph<D>::removeVertex(VertexType* pVertex)
{
	const LinearNumKey& key = pVertex->m_key;
	std::cout << "Removing " << key << "\n";

	// Update the links of this vertex to point to the merged vertex
	for(size_t dirIdx = 0; dirIdx < NUM_DIRECTIONS; ++dirIdx)
	{
		typename VertexType::EdgeCollection& currEdges = pVertex->m_edges[dirIdx];
		for(typename VertexType::EdgeCollectionIter vertexIter = currEdges.begin(); vertexIter != currEdges.end(); vertexIter++)
		{			
			// If the link has the reverse flag set, flip the direction
			extDirection vertexToChildDir = (extDirection)dirIdx;
			extDirection expectedRemoveDir = (vertexIter->reverse) ? vertexToChildDir : !vertexToChildDir;
			bool expectedRemoveReverse = vertexIter->reverse;
			
			//printVertex(vertexIter->pVertex->m_key);
			//printf("Removing %s from %s in dir %d reverseness (%d)\n", key.c_str(), vertexIter->pVertex->m_key.c_str(), expectedRemoveDir, vertexIter->reverse);
			vertexIter->pVertex->removeEdge(pVertex, expectedRemoveDir, expectedRemoveReverse);
		}
	}
	
	// Delete this vertex from the table

	// Ensure it exists
	delete pVertex;
	pVertex = NULL;
	m_vertexTable[key] = NULL;
}

//
//
//
template<typename D>
size_t DirectedGraph<D>::getDegree(const LinearNumKey& key, extDirection dir)
{
	VertexType* pVertex = findVertex(key);
	return pVertex->numEdges(dir);
}

template<typename D>
template<class Functor>
void DirectedGraph<D>::validate(Functor dataChecker)
{
	for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		VertexType* pVertex = iter->second;
		
		// validate each edge of the vertex
		for(size_t dirIdx = 0; dirIdx < NUM_DIRECTIONS; ++dirIdx)
		{
			extDirection currDir = (extDirection)dirIdx;
			typename VertexType::EdgeCollection& currEdges = pVertex->m_edges[dirIdx];
			for(typename VertexType::EdgeCollectionIter edgeIter = currEdges.begin(); edgeIter != currEdges.end(); ++edgeIter)
			{
				// Check that this edge's partner has the same orientation as it
				VertexType* pPartner = edgeIter->pVertex;
				
				extDirection expectedPartnersDir = (edgeIter->reverse) ? currDir : !currDir;
				bool expectedParentsReverse = edgeIter->reverse;
				// get the edge
				bool found;
				typename VertexType::EdgeCollectionIter partnerEdge = pPartner->getEdge(pVertex, expectedPartnersDir, expectedParentsReverse, found);
				
				if(!found)
				{
					printf("FAILED CONSISTENCY CHECK BETWEEN %s and %s (partner edge not found)\n", pVertex->m_key.c_str(), pPartner->m_key.c_str());
					exit(1);
				}
				
				
				// ensure the reverse flags are the same
				if(edgeIter->reverse != partnerEdge->reverse)
				{
					printf("FAILED CONSISTENCY CHECK BETWEEN %s and %s (reverse flags not the same)\n", pVertex->m_key.c_str(), pPartner->m_key.c_str());
					exit(1);	
				}
				
				// check that the sequences have the correct overlap
				if(!dataChecker.check(pVertex->m_data, pPartner->m_data, currDir, edgeIter->reverse))
				{
					printf("FAILED CONSISTENCY CHECK BETWEEN %s and %s (data not consistent)\n", pVertex->m_key.c_str(), pPartner->m_key.c_str());
					exit(1);	
				}
			}
		}
	}
	
	printf("Graph validated\n");
}


template<typename D>
void DirectedGraph<D>::generateComponents(VertexType* pVertex,
		extDirection dir, size_t maxCost,
		VertexComponentVector& outComponents)
{
	// Create a vertex collection for every sequence directly adjacent to this vertex
	typename VertexType::EdgeCollection& edgeCollection = pVertex->m_edges[dir];
	
	for(typename VertexType::EdgeCollectionConstIter iter = edgeCollection.begin(); iter != edgeCollection.end(); ++iter)
	{
		printf("Generating component: \n");
		// explore down this branch until maxCost is hit, accumulating all the vertices
			
		// start the branch
		VertexComponent newComp;
		newComp.first = iter->pVertex->m_key;
		outComponents.push_back(newComp);
		
		extDirection newDir = (iter->reverse) ? !dir : dir;
		accumulateVertices(iter->pVertex, newDir, 0, maxCost,
				outComponents.back().second);
	}
}

template<typename D>
void DirectedGraph<D>::accumulateVertices(VertexType* pVertex,
		extDirection dir, size_t currCost, size_t maxCost,
		VertexCollection& accumulator)
{	
	// Add this vertex
	accumulator.insert(pVertex);

	// add the cost
	currCost += costFunctor.cost(pVertex->m_data);

	if(currCost > maxCost)
	{
		return;
	}
	else
	{		
		typename VertexType::EdgeCollection& edgeCollection = pVertex->m_edges[dir];
		for(typename VertexType::EdgeCollectionConstIter iter = edgeCollection.begin(); iter != edgeCollection.end(); ++iter)
		{	
			// recursively call for each subbranch
			extDirection newDir = (iter->reverse) ? !dir : dir;
			accumulateVertices(iter->pVertex, newDir, currCost,
					maxCost, accumulator);
		}
	}
}

template<typename D>
template<class Functor>
void DirectedGraph<D>::iterativeVisit(Functor visitor)
{
	for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		visitor.visit((*iter)->m_key, (*iter)->m_data);
	}
}

template<typename D>
template<class Functor>
void DirectedGraph<D>::outputVertexConnectivity(Functor visitor) const
{
	for(VertexTableConstIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		printf("Vertex %s\n", iter->first.c_str());
		iter->second->printEdges();
		printf("Paired contigs: ");
		visitor.visit(iter->first, iter->second->m_data);
	}
}

//
// TODO: Move this logic out of this class, its way too specific
//
template<typename D>
template<class MergerFunctor>
bool DirectedGraph<D>::mergePath(const LinearNumKey& key1, const LinearNumKey& key2, extDirection parentDir, bool removeChild, bool usableChild, MergerFunctor dataMerger)
{
	VertexType* pParent = findVertex(key1);
	VertexType* pChild = findVertex(key2);

	std::cout << "Path merging " << key1 << " into " << key2 << std::endl;

	//printVertex(key1);
	//printVertex(key2);
	
	EdgeDescription parentEdgeDesc = pParent->findUniqueEdgeInDir(key2, parentDir);
	assert(merge(pParent, pChild, parentEdgeDesc.dir, parentEdgeDesc.reverse, removeChild, usableChild, dataMerger));
	
	//printVertex(key1);
	return true;
}

template<typename D>
template<class MergerFunctor>
bool DirectedGraph<D>::mergeShortestPath(const LinearNumKey& key1,
		const LinearNumKey& key2, MergerFunctor dataMerger)
{
	// Get the shortest path between the nodes
	KeyVec path;
	ShortestPathData spd;
	dijkstra(key1, spd);
	extractShortestPath(findVertex(key1), findVertex(key2), spd, path);
	mergePath(key1, path, dataMerger);
	return false;
}

//
// Compute the single-source shortest path distance to all nodes using dijkstra's algorithm
// Note that this does not consider direction so it is unsuitable to compute the djikstra shortest path if you
// are travelling in a particular direction at all times
//
template<typename D>
void DirectedGraph<D>::dijkstra(const LinearNumKey& sourceKey,
		ShortestPathData& shortestPathData)
{
	//Timer dTimer("dijkstra");

	// initiliaze infinity to a large number
	const size_t INF = 2 << 30;

	// initialize the data
	for(VertexTableConstIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		shortestPathData.distanceMap[*iter] = INF;
		shortestPathData.visitedMap[*iter] = VC_WHITE;
		shortestPathData.previousMap[*iter] = NULL;
	}
	
	VertexType* pSourceVertex = findVertex(sourceKey);
	
	VertexType* pCurrVertex = pSourceVertex;
	shortestPathData.distanceMap[pCurrVertex] = 0;
	
	bool stop = false;
	while(!stop)
	{
		shortestPathData.visitedMap[pCurrVertex] = VC_BLACK;
		
		// update all the distances of the adjacent nodes
		for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
		{
			extDirection dir = (extDirection)dirIdx;
			for(typename VertexType::EdgeCollection::iterator eIter = pCurrVertex->m_edges[dir].begin(); eIter != pCurrVertex->m_edges[dir].end(); ++eIter)
			{
				// Get the vertex to the edge points to
				VertexType* pAdjVertex = eIter->pVertex;
				
				// Get the cost to the node
				int cost = costFunctor.cost(pCurrVertex->m_data);
				
				if(shortestPathData.distanceMap[pAdjVertex] >  shortestPathData.distanceMap[pCurrVertex] + cost)
				{
					shortestPathData.distanceMap[pAdjVertex] = shortestPathData.distanceMap[pCurrVertex] + cost;
					shortestPathData.previousMap[pAdjVertex] = pCurrVertex;
					//printf("Setting prev map %s -> %s\n", pAdjVertex->m_key.c_str(), pCurrVertex->m_key.c_str());
				}
			}
		}
		
		// select the new node
		size_t minCost = INF;
		typename std::map<VertexType*, size_t>::iterator bestIter = shortestPathData.distanceMap.end();
		
		for(typename std::map<VertexType*, size_t>::iterator dIter = shortestPathData.distanceMap.begin(); dIter != shortestPathData.distanceMap.end(); ++dIter)
		{
			if(shortestPathData.visitedMap[dIter->first] != VC_BLACK)
			{
				if(dIter->second <= minCost)
				{
					bestIter = dIter;
					minCost = dIter->second;
				}
			}
		}
		
		// check if we should terminate
		if(bestIter == shortestPathData.distanceMap.end())
		{
			stop = true;
		}
		else
		{
			pCurrVertex = bestIter->first;
		}
	}
}

template<typename D>
bool DirectedGraph<D>::findSuperpaths(const LinearNumKey& sourceKey,
		extDirection dir, const KeyConstraintMap& keyConstraints,
		FeasiblePaths& superPaths,
		int maxNumPaths, int maxCompCost, int& compCost) const
{
    VertexType* pSourceVertex = findVertex(sourceKey);

    // Early exit if there are no reachable vertices
    if(keyConstraints.empty())
    {
            return false;
    }

    // Create the initial (empty) path
    VertexPath path;

	ConstrainedDFS(pSourceVertex, dir, false, keyConstraints, path,
			superPaths, 0, maxNumPaths, maxCompCost,
			compCost);
	return compCost >= maxCompCost ? false : !superPaths.empty();
}

/** Find paths through the graph that satisfy the constraints.
 * @return false if the search exited early
 */
template<typename D>
bool DirectedGraph<D>::ConstrainedDFS(VertexType* pCurrVertex,
		extDirection dir, bool rcFlip,
		const KeyConstraintMap keyConstraints,
		VertexPath currentPath, FeasiblePaths& solutions,
		size_t currLen, int maxNumPaths,
		int maxCompCost, int& visitedCount) const
{
    // Early exit if the path limit has been reached, in this case the output is invalid and should be tossed
    if ((maxNumPaths != -1 && (int)solutions.size() > maxNumPaths)
			|| ++visitedCount >= maxCompCost)
		return false;

    // Recursively explore the subbranches until either the contraints have been violated or all the constrains have been satisfied
    typename VertexType::EdgeCollection& currEdges = pCurrVertex->m_edges[dir];
    for(typename VertexType::EdgeCollection::iterator eIter = currEdges.begin(); eIter != currEdges.end(); ++eIter)
    {
        VertexType* pNextVertex = eIter->pVertex;

        // add the node to the path
        PathNode node;
        node.key = pNextVertex->m_key;
        node.isRC = eIter->reverse;
                
        // Get the relative direction and comp for the new node
        extDirection relativeDir = EdgeDescription::getRelativeDir(dir, eIter->reverse);
        bool relativeRC = rcFlip ^ eIter->reverse;

        VertexPath newPath = currentPath;
        newPath.push_back(node);

        // Update the constraints set
        KeyConstraintMap newConstraints = keyConstraints;

        // Make a copy of the constraints
        typename KeyConstraintMap::iterator constraintIter = newConstraints.find(pNextVertex->m_key);
        if(constraintIter != newConstraints.end())
        {
        		if(constraintIter->second.isConstraintMet(currLen, relativeRC))
        		{
	                // This vertex is a valid constraint, remove it from the set
	                newConstraints.erase(constraintIter);
        		}
        		/*else
        		{
        			constraintIter->second.printStatus(currLen, relativeRC);
        		}*/
        }

        // Have all the constraints been satisfied?
        if(newConstraints.empty())
        {
                // this path is valid
                solutions.push_back(newPath);
                continue;
        }
        else
        {
        	
	            // update the path length
	            size_t newLength = currLen + costFunctor.cost(pNextVertex->m_data);    
	            
                // Check if this path can possibly support all the constraints
                bool constraintViolated = false;
                for(typename KeyConstraintMap::iterator cIter = newConstraints.begin(); cIter != newConstraints.end(); ++cIter)
                {
                        //std::cout << "a " << cIter->second << " , " << (int)newLength << " result " << (cIter->second < (int)newLength) << "\n";
                	 	//if(cIter->second.distance < (int)newLength)
                		if(!cIter->second.isConstraintPossible(newLength))
                        {
                                //this constraint is violated, the branch will be explored no further
                                constraintViolated = true;
                                break;
                        }
                }

                if(constraintViolated)
                {
                       continue;
                }

                //printf("path length: %zu, num contraints: %zu\n", newLength, newConstraints.size());

                // recurse
				if (!ConstrainedDFS(pNextVertex, relativeDir,
						relativeRC,
						newConstraints, newPath, solutions,
						newLength, maxNumPaths,
						maxCompCost, visitedCount))
					return false;
        }
    }
	return true;
}

//
// Return the minimum possible path length that will contain every vertex in the set
//
template<typename D>
size_t  DirectedGraph<D>::getMinPathLength(
		const VertexPtrSet& vertexSet)
{
	// The minimum possible path length has the longest node as the terminal
	
	// Sum the overlaps for every node to get the total path length
	size_t pathLength = 0;
	size_t maxCost = 0;
	for(typename VertexPtrSet::iterator iter = vertexSet.begin(); iter != vertexSet.end(); ++iter)
	{
		// add the cost of going through this vertex
		size_t vertexCost = costFunctor.cost((*iter)->m_data);
		pathLength += vertexCost;
		
		if(vertexCost > maxCost)
		{
			maxCost = vertexCost;
		}
	}
	
	// Subtract the largest cost
	pathLength -= maxCost;
	return pathLength;
}

//
//
//
template<typename D>
void DirectedGraph<D>::extractShortestPath(VertexType* pSource, VertexType* pTarget, ShortestPathData& shortestPathData, KeyVec& path)
{
	VertexType* pCurrVertex = pTarget;
	while(pCurrVertex != pSource)
	{
		path.push_back(pCurrVertex->m_key);
		pCurrVertex = shortestPathData.previousMap[pCurrVertex];
	}
	
	// Reverse the vector
	std::reverse(path.begin(), path.end());
}

template<typename D>
size_t DirectedGraph<D>::calculatePathLength(const VertexPath& path)
	const
{
	size_t len = 0;
	for(typename VertexPath::const_iterator iter = path.begin(); iter != path.end() - 1; ++iter)
	{
		len += costFunctor.cost(getDataForVertex(iter->key));
	}
	return len;
}

template<typename D>
void DirectedGraph<D>::makeDistanceMap(const VertexPath& path,
		std::map<LinearNumKey, int>& distanceMap) const
{
	// the path distance to a node is the distance that walks through all the nodes leading to it
	// the first node in a path therefore has a distance of 0 by def
	size_t distance = 0;
	for(typename VertexPath::const_iterator iter = path.begin(); iter != path.end(); ++iter)
	{
		distanceMap[iter->key] = distance;
		int currCost = costFunctor.cost(getDataForVertex(iter->key));
		distance += currCost;		
	}
	return;	
}
