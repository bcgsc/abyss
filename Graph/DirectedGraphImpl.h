//
//
// Vertex Implemenation
//
//

//
//
//
template<typename K, typename D>
size_t Vertex<K,D>::findEdgeIndex(const K& key)
{
	for(size_t idx = 0; idx <= 1; ++idx)
	{
		const EdgeCollection& currEdgeSet = m_edges[idx];
		for(EdgeCollectionConstIter edgeIter = currEdgeSet.begin(); edgeIter != currEdgeSet.end(); ++edgeIter)
		{
			if(edgeIter->first->m_key == key)
			{
				return idx;
			}
		}
	}
	
	// not found code
	return 2;
}	

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
		std::cout << "Vertex data: " << pVertex->m_data.seq << std::endl;
	}
}

//
// Attempt to reduce the data set using paired reads
//
template<typename K, typename D>
template<class DataCostFunctor, class ResolveFunctor, class DataMerger>
size_t DirectedGraph<K,D>::reducePaired(size_t maxCost, DataCostFunctor& dataCost, ResolveFunctor& resolver, DataMerger& merger)
{
	size_t numMerged = 0;
	
	//while(!stop)

	for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
	{
		// hack
		if(iter->second->m_data.seq.length() > 200)
		{
			for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
			{
				while(attemptResolve(iter->first, (extDirection)dirIdx, maxCost, dataCost, resolver, merger));
			}
		}
	}

	removeTransitivity(merger);
	return numMerged;
		
/*
	size_t numMerged = 0;
	
	bool stop = false;
	while(!stop)
	{
MERGESTART:
		bool merged = false;
		for(VertexTableIter iter = m_vertexTable.begin(); iter != m_vertexTable.end(); ++iter)
		{
			for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
			{
				merged = attemptResolve(iter->first, (extDirection)dirIdx, maxCost, dataCost, resolver, merger);
				if(merged)
				{
					removeTransitivity(merger);
					printf("num nodes now: %zu\n", m_vertexTable.size());
					goto MERGESTART;
				}
			}
		}
		
		if(!merged)
		{
			stop = true;
		}
	}
	return numMerged;
*/
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
				bool reverse = currEdges.begin()->reverse;
				

									
				// attempt the merge
				if(merge(pVertex, pPartner, (extDirection)idx, reverse, dataMerger))
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


				
	
	// Find the link to the child
	for(size_t idx = 0; idx < NUM_DIRECTIONS; ++idx)
	{
		bool found;
		
		for(size_t revIdx = 0; revIdx <= 1; ++revIdx)
		{
			bool reverse = (revIdx == 1);
			pParent->getEdge(pChild, (extDirection)idx, reverse, found);
			
			if(found)
			{
				return merge(pParent, pChild, (extDirection)idx, reverse, dataMerger);
			}	
		}
	}
	return false;
}

template<typename K, typename D>
template<class Functor>
bool DirectedGraph<K,D>::merge(VertexType* pParent, VertexType* pChild, const extDirection parentsDir, const bool parentsReverse, Functor dataMerger)
{
	// Determine the relative orientation of the nodes
	// If they point towards each other, the child vertex must be flipped around
	// ---SENSE---> <---SENSE = FLIP
	// ----SENSE---> <---ANTISENSE = OK
	
	/*
	printf("trying to merge %s to %s (%d %d)\n", pParent->m_key.c_str(), pChild->m_key.c_str(), parentsDir, parentsReverse);
	printf("pre merge\n");
	printVertex(pParent->m_key);
	printVertex(pChild->m_key);
	*/
	
	K parentKey = pParent->m_key;
	K childKey = pChild->m_key;
	
	// check if this node has been merged before
	if(pParent->m_mergeRecord.find(childKey) != pParent->m_mergeRecord.end())
	{
		//do not merge, loop
		return false;
	}
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
	dataMerger.merge(pParent->m_data, pChild->m_data, (extDirection)parentsDir, parentsReverse);
	
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
	
	// if the child now has no links in this direction, it can be removed as it is redundant
	if(pChild->m_edges[expectedChildsDir].empty())
	{
		// remove the vertex and update all the links to point to the merged node
		removeVertex(pChild);
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
void DirectedGraph<K,D>::removeVertex(VertexType* pVertex)
{
	const K& key = pVertex->m_key;
	//printf("Removing %s\n", key.c_str());

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
	
	//printf("deleting %s\n", key.c_str());

	
	delete pVertex;
	pVertex = NULL;
	
	m_vertexTable.erase(iter);
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
template<class DataCostFunctor, class ResolveFunctor, class DataMerger>
bool DirectedGraph<K,D>::attemptResolve(const K& key, extDirection dir, size_t maxCost, DataCostFunctor& dataCost, ResolveFunctor& resolver, DataMerger& merger)
{
	// Find the vertex to resolve
	VertexType* pVertex = findVertex(key);
	
	printf("Before merge\n");
	printVertex(key);
	
	// Generate the vertex collections for the adjacent sequences
	VertexComponentVector extensionComponents;
	generateComponents(pVertex, dir, maxCost, extensionComponents, dataCost);
	
	// Generate the data components which the resolver will operate on
	DataComponents dataComponents;

	size_t count = 0;
	for(typename VertexComponentVector::iterator compIter = extensionComponents.begin(); compIter != extensionComponents.end(); ++compIter)
	{
		DataCollection dataCollection;
		printf("Component %zu (%zu)\n", count++, compIter->second.size());
		for(typename VertexCollection::iterator vertIter = compIter->second.begin(); vertIter != compIter->second.end(); ++vertIter)
		{
			printf("	key %s len %zu\n", (*vertIter)->m_key.c_str(), (*vertIter)->m_data.seq.length());
			dataCollection.push_back((*vertIter)->m_data);
		}
		
		dataComponents.push_back(dataCollection);
	}
	
	// The resolver will attempt to choose a component index to merge to
	int index = resolver.resolve(pVertex->m_data, dir, dataComponents);
	
	if(index == -1)
	{
		// could not be resolved
		return false;
	}
	else
	{
		// get the key of the supported branch
		const K selectedKey = extensionComponents[index].first;
		printf("Selected key for merge: %s\n", selectedKey.c_str());
		
		// perform the merge
		
		// TODO: FIX THIS
		return mergeWrapper(pVertex->m_key, selectedKey, merger);
	}
	
}

template<typename K, typename D>
template<class DataCostFunctor>
void DirectedGraph<K,D>::generateComponents(VertexType* pVertex, extDirection dir, size_t maxCost, VertexComponentVector& outComponents, DataCostFunctor& dataCost)
{
	// Create a vertex collection for every sequence directly adjacent to this vertex
	typename VertexType::EdgeCollection& edgeCollection = pVertex->m_edges[dir];
	for(typename VertexType::EdgeCollectionConstIter iter = edgeCollection.begin(); iter != edgeCollection.end(); ++iter)
	{
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


/*
template<typename K, typename D>
template<class Functor>
void DirectedGraph<T,D>::bfsVisit(Functor f)
{
	T def("ATAGAC");
	f(def);
}*/
