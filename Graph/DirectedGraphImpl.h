//
//
// Vertex Implemenation
//
//

//
//
//
template<typename K, typename D>
void Vertex<K,D>::addEdge(VertexType* pNode, extDirection dir, bool reverse) 
{
	EdgeData data;
	data.pVertex = pNode;
	data.reverse = reverse;
	m_edges[dir].insert(data);
}

//
//
//
template<typename K, typename D>
void Vertex<K,D>::removeEdge(VertexType* pNode, extDirection dir, bool reverse)
{
	// Remove the node
	EdgeData data;
	data.pVertex = pNode;
	data.reverse = reverse;
		
	EdgeCollection& currEdgeSet = m_edges[dir];
	EdgeCollectionIter iter = currEdgeSet.find(data);
	assert(iter != currEdgeSet.end());
	currEdgeSet.erase(iter);
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
	
	//printf("%s looking for %s in dir %d with rev %d (num links %zu)\n", m_key.c_str(), pNode->m_key.c_str(), dir, reverse, m_edges[dir].size());
	
	EdgeCollectionIter iter = m_edges[dir].find(data);
	if(iter != m_edges[dir].end())
	{
		found = true;
	}
	else
	{
		found = false;
	}
	return iter;
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
				if(!found)
				{
					ret.dir = (extDirection)dirIdx;
					ret.reverse = edgeIter->reverse;
					found = true;
				}
				else
				{
					printf("Error: Tried to get an edge without checking for uniqueness!\n");
					assert(false);
				}
			}
		}
	}

	if(!found)
	{
		printf("Could not find edge %s for vertex %s\n", key.c_str(), m_key.c_str());
		assert(found);
	}
	return ret;
}

//
//
//
template<typename K, typename D>
EdgeDescription Vertex<K,D>::findUniqueEdgeInDir(const K& key, extDirection dir)
{
	bool found = false;
	EdgeDescription ret;

	const EdgeCollection& currEdgeSet = m_edges[dir];
	for(EdgeCollectionConstIter edgeIter = currEdgeSet.begin(); edgeIter != currEdgeSet.end(); ++edgeIter)
	{
		if(edgeIter->pVertex->m_key == key)
		{
			if(!found)
			{
				ret.dir = dir;
				ret.reverse = edgeIter->reverse;
				found = true;
			}
			/*else
			{
				printf("Error: Tried to get an edge without checking for uniqueness!\n");
				assert(false);
			}*/
		}
	}

	if(!found)
	{
		printf("Could not find edge %s for vertex %s\n", key.c_str(), m_key.c_str());
		assert(found);
	}
	
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
template<typename K, typename D>
DirectedGraph<K,D>::~DirectedGraph()
{
	// delete all the vertices in the table
	for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		delete iter->second;	
	}
}

template<typename K, typename D>
void DirectedGraph<K,D>::addEdge(const K& parent, const K& child, extDirection dir, bool reverse)
{	
	// Get iterators to the nodes
	// If they do not exist, get will add them
	VertexType* pParentVertex = findVertex(parent);
	assert(pParentVertex != NULL);
	
	VertexType* pChildVertex = findVertex(child);
	assert(pChildVertex != NULL);
	
	// Add the link
	pParentVertex->addEdge(pChildVertex, dir, reverse);
}

template<typename K, typename D>
typename DirectedGraph<K,D>::VertexType* DirectedGraph<K,D>::addVertex(const K& key, const D& data)
{
	//assert(findVertex(key) == NULL);
	
	// Create new vertex
	VertexType* ptr = new VertexType(key, data);
	
	// Insert into the table
	m_vertexTable[key] = ptr;
	
	return ptr;
}

template<typename K, typename D>
typename DirectedGraph<K,D>::VertexType* DirectedGraph<K,D>::findVertex(const K& key) const
{
	// Look up the key in the vertex table
	VertexTableConstIter iter = m_vertexTable.find(key);
	
	if(iter != m_vertexTable.end())
	{
		return iter->second;
	}
	else
	{
		return NULL;
	}
}

// Count all the edges in all the nodes
template<typename K, typename D>
size_t DirectedGraph<K,D>::countEdges() const
{
	size_t sum = 0;
	for(VertexTableConstIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		sum += iter->second->countEdges();
	}
	return sum;
}

// Print out a node
template<typename K, typename D>
void DirectedGraph<K,D>::printVertex(const K& key, bool printData) const
{
	// find the vertex
	VertexType* pVertex = findVertex(key);
	
	// should not be null
	assert(pVertex != NULL);
	
	// Print the vertex's key
	std::cout << "Vertex key: " << pVertex->m_key << std::endl;
	pVertex->printEdges();
	
	if(printData)
	{
		// debug, not template friendly
		std::cout << "Vertex data: " << pVertex->m_data.m_seq << std::endl;
	}
}

//
// Attempt to reduce the data set using paired reads
//
template<typename K, typename D>
template<class ResolveFunctor>
size_t DirectedGraph<K,D>::reducePaired(ResolveFunctor& resolver)
{
	size_t numMerged = 0;
	
	for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		if(iter->second->m_data.m_seq.length() > 10000)
		{
			while(resolver.resolve(this, iter->first));
		}
	}

	return numMerged;
		
}

//
// Remove transitivity in the data set by merging and 
//
template<typename K, typename D>
template<class Functor>
size_t DirectedGraph<K,D>::removeTransitivity(Functor dataMerger)
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
				
				// Can merge edge
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
template<typename K, typename D>
template<class Functor>
bool DirectedGraph<K,D>::mergeWrapper(const K& key1, const K& key2, Functor dataMerger)
{
	VertexType* pParent = findVertex(key1);
	VertexType* pChild = findVertex(key2);

	extDirection dir;
	bool reverse;
	pParent->findUniqueEdge(key2, dir, reverse);
	
	extDirection childDir = EdgeDescription::getTwinDir(dir, reverse);
	
	// Should the child be removed?
	bool removeChild = false;
	
	// is there only one edge between the parent and child?
	if(pChild->numEdges(childDir) == 1)
	{
		removeChild = true;
	}
	
	return merge(pParent, pChild, dir, reverse, removeChild, dataMerger);
}

//
// Append a copy of the child vertex into the parent and update all the links accordingly
//
template<typename K, typename D>
template<class Functor>
bool DirectedGraph<K,D>::merge(VertexType* pParent, VertexType* pChild, const extDirection parentsDir, const bool parentsReverse, bool removeChild, Functor dataMerger)
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
	
	K parentKey = pParent->m_key;
	K childKey = pChild->m_key;
	
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
	dataMerger.merge(pParent->m_key, pParent->m_data, pChild->m_key, pChild->m_data, (extDirection)parentsDir, parentsReverse, removeChild);
	
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
		removeVertex(pChild, dataMerger);
	}
	
	/*
	printf("post merge\n");
	printVertex(pParent->m_key);
	//printVertex(pPartner->m_key);	
	*/
	
	pParent->m_mergeRecord.insert(childKey);
	
	// merged, return true
	return true;
}

template<typename K, typename D>
template<class DataFunctor>
void DirectedGraph<K,D>::removeVertex(VertexType* pVertex, DataFunctor functor)
{
	const K& key = pVertex->m_key;
	printf("Removing %s\n", key.c_str());

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
	VertexTableIter iter = m_vertexTable.find(key);
	assert(iter != m_vertexTable.end());
	
	// Call the functor to indicate this vertex is being deleted
	functor.deleteCallback(key, pVertex->m_data);
	//printf("deleting %s\n", key.c_str());

	
	delete pVertex;
	pVertex = NULL;
	
	m_vertexTable.erase(iter);
}

//
//
//
template<typename K, typename D>
size_t DirectedGraph<K,D>::getDegree(const K& key, extDirection dir)
{
	VertexType* pVertex = findVertex(key);
	return pVertex->numEdges(dir);
}

template<typename K, typename D>
template<class Functor>
void DirectedGraph<K,D>::validate(Functor dataChecker)
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


template<typename K, typename D>
template<class DataCostFunctor>
void DirectedGraph<K,D>::generateComponents(VertexType* pVertex, extDirection dir, size_t maxCost, VertexComponentVector& outComponents, DataCostFunctor& dataCost)
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
		accumulateVertices(iter->pVertex, newDir, 0, maxCost, outComponents.back().second, dataCost);
	}
}

template<typename K, typename D>
EdgeDescription DirectedGraph<K,D>::getUniqueEdgeDesc(const K& key1, const K& key2, extDirection parentDir)
{
    VertexType* pSourceVertex = findVertex(key1);
    return pSourceVertex->findUniqueEdgeInDir(key2, parentDir);
}

template<typename K, typename D>
template<class DataCostFunctor>
void DirectedGraph<K,D>::accumulateVertices(VertexType* pVertex, extDirection dir, size_t currCost, size_t maxCost, VertexCollection& accumulator, DataCostFunctor& dataCost)
{	
	// Add this vertex
	accumulator.insert(pVertex);
	
	// add the cost
	currCost += dataCost.cost(pVertex->m_data);
	
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
			accumulateVertices(iter->pVertex, newDir, currCost, maxCost, accumulator, dataCost);
		}
	}
}

template<typename K, typename D>
template<class Functor>
void DirectedGraph<K,D>::iterativeVisit(Functor visitor)
{
	for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		visitor.visit(iter->first, iter->second->m_data);
	}
}

template<typename K, typename D>
template<class Functor>
void DirectedGraph<K,D>::outputVertexConnectivity(Functor visitor) const
{
	for(VertexTableConstIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		printf("Vertex %s\n", iter->first.c_str());
		iter->second->printEdges();
		printf("Paired contigs: ");
		visitor.visit(iter->first, iter->second->m_data);
	}
}

template<typename K, typename D>
template<class MergerFunctor>
bool DirectedGraph<K,D>::mergePath(const K& key1, const K& key2, extDirection parentDir, bool removeChild, MergerFunctor dataMerger)
{
	VertexType* pParent = findVertex(key1);
	VertexType* pChild = findVertex(key2);

	printf("Path merging %s with %s\n", key1.c_str(), key2.c_str());
	printVertex(key1);
	printVertex(key2);
	
	EdgeDescription parentEdgeDesc = pParent->findUniqueEdgeInDir(key2, parentDir);
	assert(merge(pParent, pChild, parentEdgeDesc.dir, parentEdgeDesc.reverse, removeChild, dataMerger));
	
	printVertex(key1);
	return true;
}

template<typename K, typename D>
template<class DataCostFunctor, class MergerFunctor>
bool DirectedGraph<K,D>::mergeShortestPath(const K& key1, const K& key2, DataCostFunctor costFunctor, MergerFunctor dataMerger)
{
	// Get the shortest path between the nodes
	KeyVec path;
	ShortestPathData spd;
	dijkstra(key1, spd, costFunctor);
	extractShortestPath(findVertex(key1), findVertex(key2), spd, path);
	mergePath(key1, path, dataMerger);
	return false;
}

//
// Generate the shortest directed path to any node in the set by using a greedy algorithm
//
template<typename K, typename D>
template<class DataCostFunctor>
void DirectedGraph<K,D>::greedyDirectedPath(const K& sourceKey, extDirection dir, KeySet& /*terminals*/, ShortestPathData& shortestPathData, DataCostFunctor& /*costFunctor*/)
{
	Timer t("BFS");
	// initiliaze infinity to a large number
	const size_t INF = 2 << 30;

	// initialize the data
	for(VertexTableConstIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		shortestPathData.distanceMap[iter->second] = INF;
		shortestPathData.visitedMap[iter->second] = VC_WHITE;
		shortestPathData.previousMap[iter->second] = NULL;
	}
	
	VertexType* pSourceVertex = findVertex(sourceKey);
	
	shortestPathData.distanceMap[pSourceVertex] = 0;
	shortestPathData.visitedMap[pSourceVertex] = VC_GRAY;
	
	std::queue<VertexDirPair> visitQueue;
	VertexDirPair sourcePair;
	sourcePair.pVertex = pSourceVertex;
	sourcePair.dir = dir;
	
	visitQueue.push(sourcePair);
	
	while(!visitQueue.empty())
	{
		VertexDirPair currentPair = visitQueue.front();
		VertexType* pCurrVertex = currentPair.pVertex;
		extDirection currentDir = currentPair.dir;
		
		visitQueue.pop();
		
		typename VertexType::EdgeCollection& currEdges = pCurrVertex->m_edges[currentDir];
		for(typename VertexType::EdgeCollection::iterator eIter = currEdges.begin(); eIter != currEdges.end(); ++eIter)
		{
			if(shortestPathData.visitedMap[eIter->pVertex] == VC_WHITE)
			{
				VertexType* pAdjVertex = eIter->pVertex;
				shortestPathData.visitedMap[pAdjVertex] = VC_GRAY;
				shortestPathData.distanceMap[pAdjVertex] = shortestPathData.distanceMap[pCurrVertex] + 1;
				shortestPathData.previousMap[pAdjVertex] = pCurrVertex;
				
				// Get the relative direction
				extDirection relativeDir = EdgeDescription::getRelativeDir(currentDir, eIter->reverse);
				VertexDirPair np;
				np.pVertex = pAdjVertex;
				np.dir = relativeDir;
				
				visitQueue.push(np);
			}
		}
		
		shortestPathData.visitedMap[pCurrVertex] = VC_BLACK;
		shortestPathData.directionMap[pCurrVertex] = currentDir;
	}
}

//
// Compute the single-source shortest path distance to all nodes using dijkstra's algorithm
// Note that this does not consider direction so it is unsuitable to compute the djikstra shortest path if you
// are travelling in a particular direction at all times
//
template<typename K, typename D>
template<class DataCostFunctor>
void DirectedGraph<K,D>::dijkstra(const K& sourceKey, ShortestPathData& shortestPathData, DataCostFunctor& costFunctor)
{
	//Timer dTimer("dijkstra");

	// initiliaze infinity to a large number
	const size_t INF = 2 << 30;

	// initialize the data
	for(VertexTableConstIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		shortestPathData.distanceMap[iter->second] = INF;
		shortestPathData.visitedMap[iter->second] = VC_WHITE;
		shortestPathData.previousMap[iter->second] = NULL;
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

template<typename K, typename D>
template<class DataCostFunctor>
void DirectedGraph<K,D>::computeMinimalPath(const K& sourceKey, extDirection dir, VertexPtrSet vertexConstraints, const size_t maxPathLen, VertexPath& outpath, DataCostFunctor& costFunctor)
{
	(void)maxPathLen;
	VertexType* pSourceVertex = findVertex(sourceKey);
	VertexType* pCurrVertex = pSourceVertex;	
	
	// Create the initial (empty) path
	VertexPath path;
	
	size_t minLen = getMinPathLength(vertexConstraints, costFunctor);
	printf("Min path length: %zu\n", minLen);
	
	// Recursively search each sub branch for a feasible solution
	FeasiblePaths solutions;
	ConstrainedDFS(pCurrVertex, dir, vertexConstraints, path, solutions, 0, maxPathLen, costFunctor);
	
	printf("%zu feasible paths found\n", solutions.size());	
	
	size_t bestLen = maxPathLen;
	typename FeasiblePaths::const_iterator bestIter = solutions.end();
	
	// Select a path with a minimum length (TODO: a better heuristic would be a max likelihood that the path is correct)
	for(typename FeasiblePaths::const_iterator iter = solutions.begin(); iter != solutions.end(); ++iter)
	{
		printVertexPath(*iter, costFunctor);	
		size_t pathLen = calculatePathLength(*iter, costFunctor);
		
		if(pathLen < bestLen)
		{
			bestLen = pathLen;
			bestIter = iter;
		}
	}
	
	// Was a best path chosen?
	if(bestIter != solutions.end())
	{
		outpath = *bestIter;
	}
}

template<typename K, typename D>
template<class DataCostFunctor>
void DirectedGraph<K,D>::ConstrainedDFS(VertexType* pCurrVertex, extDirection dir, VertexPtrSet vertexConstraints, 
										VertexPath currentPath, FeasiblePaths& solutions,
										size_t currLen, const size_t maxPathLen, DataCostFunctor& costFunctor)
{
	// Recursively explore the subbranches until either the contraints have been violated or all the constrains have been satisfied
	typename VertexType::EdgeCollection& currEdges = pCurrVertex->m_edges[dir];
	for(typename VertexType::EdgeCollection::iterator eIter = currEdges.begin(); eIter != currEdges.end(); ++eIter)
	{
		VertexType* pNextVertex = eIter->pVertex;
		
		// add the node to the path
		VertexPath newPath = currentPath;
		newPath.push_back(pNextVertex);

		// Update the constraints set
		VertexPtrSet newConstraints = vertexConstraints;
		newConstraints.erase(pNextVertex);
		
		// Have all the constraints been satisfied?
		if(newConstraints.empty())
		{
			// this path is valid
			solutions.push_back(newPath);
			return;
		}
		else
		{
			// update the path length
			size_t newLength = currLen + costFunctor.cost(pNextVertex->m_data);
			
			// get the minimum path length of the remaining nodes in the constraint set
			size_t minLenRemaining = getMinPathLength(newConstraints, costFunctor);
			
			//printf("path length: %zu, minRemaining: %zu num contraints: %zu\n", newLength, minLenRemaining, newConstraints.size());
			
			// check if the branch is no longer feasible
			if(newLength + minLenRemaining > maxPathLen)
			{
				// the solution is not feasible, explore this branch no further
				continue;
			}
			else
			{
				// Get the relative direction for the new node
				extDirection relativeDir = EdgeDescription::getRelativeDir(dir, eIter->reverse);	
				
				// recurse
				ConstrainedDFS(pNextVertex, relativeDir, newConstraints, newPath, solutions, newLength, maxPathLen, costFunctor);
			}
		}
	}
}

//
// Return the minimum possible path length that will contain every vertex in the set
//
template<typename K, typename D>
template<class DataCostFunctor>
size_t  DirectedGraph<K,D>::getMinPathLength(const VertexPtrSet& vertexSet, DataCostFunctor costFunctor)
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
// Find a superpath that includes all the nodes in the reachable set
//
template<typename K, typename D>
template<class DataCostFunctor>
bool DirectedGraph<K,D>::findSuperpath(const K& sourceKey, extDirection dir, const KeySet& reachableSet, KeyVec& superPath, DataCostFunctor& costFunctor)
{
	(void)dir;
	(void)superPath;
		
	VertexType* pSourceVertex = findVertex(sourceKey);
	VertexType* pCurrVertex = pSourceVertex;
	KeySet visitSet = reachableSet;
	extDirection currentDir = dir;

	// Early exit if there are no reachable vertices
	if(reachableSet.empty())
	{
		return false;
	}
	
	// Convert the key set to a vertex set
	VertexPtrSet constraints;
	
	for(typename KeySet::iterator iter = reachableSet.begin(); iter != reachableSet.end(); ++iter)
	{
		constraints.insert(findVertex(*iter));
	}
	
	// Compute the minimial path that reaches all the vertices within the specified cost
	VertexPath chosenPath;
	computeMinimalPath(pCurrVertex->m_key, currentDir, constraints, 250, chosenPath, costFunctor);
	
	if(!chosenPath.empty())
	{
		printf("Best Path:\n");
		printVertexPath(chosenPath, costFunctor);
		
		// Convert the path to a keyvec
		for(typename VertexPath::const_iterator pathIter = chosenPath.begin(); pathIter != chosenPath.end(); ++pathIter)
		{
			superPath.push_back((*pathIter)->m_key);
		}
		return true;
	}
	else
	{
		printf("No path found\n");
		return false;
	}

	
}

//
//
//
template<typename K, typename D>
void DirectedGraph<K,D>::extractShortestPath(VertexType* pSource, VertexType* pTarget, ShortestPathData& shortestPathData, KeyVec& path)
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

template<typename K, typename D>
template<class DataCostFunctor>
void DirectedGraph<K,D>::printVertexPath(const VertexPath& path, DataCostFunctor costFunctor)
{
	size_t length = calculatePathLength(path, costFunctor);
	
	
	printf("PathLen: %zu [ ", length);
	for(typename VertexPath::const_iterator iter = path.begin(); iter != path.end(); ++iter)
	{
		printf("%s ", (*iter)->m_key.c_str());
	}
	printf("] \n");
}

template<typename K, typename D>
template<class DataCostFunctor>
size_t DirectedGraph<K,D>::calculatePathLength(const VertexPath& path, DataCostFunctor costFunctor)
{
	size_t len = 0;
	for(typename VertexPath::const_iterator iter = path.begin(); iter != path.end() - 1; ++iter)
	{
		len += costFunctor.cost((*iter)->m_data);
	}
	return len;
}

/*
template<typename K, typename D>
template<class Functor>
void DirectedGraph<T,D>::bfsVisit(Functor f)
{
	T def("ATAGAC");
	f(def);
}*/
