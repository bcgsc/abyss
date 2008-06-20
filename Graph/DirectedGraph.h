#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H

#include "Timer.h"
#include "SetOperations.h"
#include <queue>

enum VisitColor
{
	VC_WHITE,
	VC_GRAY,
	VC_BLACK
};

struct EdgeDescription
{
	extDirection dir;
	bool reverse;
	
	// Get the twin direction of the specified edge
	// If the nodes are the same comp (reverse == false) then it is the opposite of dir, otherwise it is die
	static extDirection getTwinDir(extDirection refDir, bool reverse) { return ((reverse) ? refDir : !refDir); }
	
	// Get the relative direction of the specified edge
	// If the node has the same comp, it is in the same dir
	static extDirection getRelativeDir(extDirection refDir, bool reverse) { return ((!reverse) ? refDir : !refDir); }
	
	EdgeDescription makeRelativeDesc() 
	{ 
		EdgeDescription ret;
		ret.reverse = this->reverse;
		ret.dir = getRelativeDir(this->dir, this->reverse); 
		return ret;
	}  
};

template<typename K, typename D>
struct Vertex
{
	Vertex(const K& k, const D& d) : m_key(k), m_data(d) { }
	K m_key;
	D m_data;
		
	// typedefs for convienance
	typedef Vertex<K,D> VertexType;
	typedef std::set<K> MergeRecord;
	
	struct EdgeData
	{
		VertexType* pVertex;
		bool reverse;
		
		bool operator<(const EdgeData& e2) const 
		{
			//int cmp = this->pVertex - e2.pVertex;
			int cmp = strcmp(this->pVertex->m_key.c_str(), e2.pVertex->m_key.c_str());
			
			if (cmp == 0)
			{
				return this->reverse < e2.reverse;
			}
			return (cmp < 0);	  	
		} 
	  
	};
		
	typedef typename std::set<EdgeData> EdgeCollection;
	typedef typename EdgeCollection::iterator EdgeCollectionIter;
	typedef typename EdgeCollection::const_iterator EdgeCollectionConstIter;
	
	
	
	// Compare operators (needed?)
	inline bool operator==(const Vertex& other) const { return m_data == other.m_data; }
	inline bool operator<(const Vertex& other) const { return m_data < other.m_data; }
	
	// count the edges the vertex contains
	size_t countEdges() const { return (numEdges(SENSE) + numEdges(ANTISENSE)); }
	
	// get the number of edges in a certain direction
	size_t numEdges(extDirection dir ) const { return m_edges[dir].size(); }
		
	// add an edge to the vertex in the specified direction
	void addEdge(VertexType* pNode, extDirection dir, bool reverse);
	
	// remove an edge
	void removeEdge(VertexType* pNode, extDirection dir, bool reverse);
	
	// find an edge, will assert if the vertex has more than one edge pointing to the key
	EdgeDescription findUniqueEdge(const K& key);
	
	// find an edge, will assert if the vertex has more than one edge pointing to the key
	EdgeDescription findUniqueEdgeInDir(const K& key, extDirection dir);	
	
	// check if the edge with the specified directionality exists
	bool edgeExists(const K& key, extDirection dir, bool reverse);
	
	// check if the described edge is unique
	bool isEdgeUnique(VertexType* pNode, extDirection dir, bool reverse);
	
	// get an edge
	EdgeCollectionIter getEdge(VertexType* pNode, extDirection dir, bool reverse, bool& found);
	
	// detect simple cycle
	bool detectSimpleCycle();

	// Print the links of the graph
	void printLinks(const EdgeCollection& collection) const;
	
	// Print the edges of this node
	void printEdges() const;
		
	EdgeCollection m_edges[NUM_DIRECTIONS];
	MergeRecord m_mergeRecord;
	
};

template<typename K, typename D>
class DirectedGraph 
{
	// Lots of typedefs
	typedef typename std::map<K, Vertex<K, D>* > VertexTable;
	typedef typename VertexTable::iterator VertexTableIter;
	typedef typename VertexTable::const_iterator VertexTableConstIter;
	typedef Vertex<K, D> VertexType;
	
	typedef std::set<VertexType*> VertexCollection;
	typedef std::pair<K, VertexCollection> VertexComponent;
	typedef std::vector<VertexComponent> VertexComponentVector;
	
	typedef std::set<K> KeySet;
	typedef std::vector<K> KeyVec;
	typedef std::set<VertexType*> VertexPtrSet;
	
	struct VertexDirPair
	{
		VertexType* pVertex;
		extDirection dir;
	};
	
	typedef std::vector<VertexType*> VertexPath;
	typedef std::vector<VertexPath> FeasiblePaths;
	
	typedef std::vector<D*> DataCollection;
	typedef std::vector<DataCollection> DataComponents;
	
	typedef std::map<VertexType*, size_t> DistanceMap;
	typedef std::map<VertexType*, VisitColor> VisitedMap;
	typedef std::map<VertexType*, VertexType*> PreviousMap;
	typedef std::map<VertexType*, extDirection> DirectionMap;
	
	struct ShortestPathData
	{
		DistanceMap distanceMap;
		VisitedMap visitedMap;
		PreviousMap previousMap;
		DirectionMap directionMap;
	};
	
	public:
		DirectedGraph() { };
		~DirectedGraph();
		
		// get the data pointer
		const D& getDataForVertex(const K& key) const { VertexType* pVertex = findVertex(key); assert(pVertex != NULL); return pVertex->m_data; }
		
		// add edge
		void addEdge(const K& parent, const K& child, extDirection dir, bool reverse);
		
		// add vertex
		VertexType* addVertex(const K& key, const D& data);
		
		// remove vertex
		void removeVertex(VertexType* pVertex);
		
		// reduce the graph with paired data		
		template<class ResolveFunctor>
		size_t reducePaired(ResolveFunctor& resolver);
		
		// attempt to resolve 
		template<class DataCostFunctor, class ResolveFunctor, class DataMerger>
		bool attemptResolve(const K& key, extDirection dir, size_t maxCost, DataCostFunctor& dataCost, ResolveFunctor& resolver, DataMerger& merger);
		
		template<class DataCostFunctor>
		void generateComponents(VertexType* pVertex, extDirection dir, size_t maxCost, VertexComponentVector& outComponents, DataCostFunctor& dataCost);

		template<class DataCostFunctor>
		void accumulateVertices(VertexType* pVertex, extDirection dir, size_t currCost, size_t maxCost, VertexCollection& accumulator, DataCostFunctor& dataCost);
		
		template<class DataCostFunctor>
		bool findSuperpaths(const K& sourceKey, extDirection dir, const KeySet& reachableSet, std::vector<KeyVec>& solutions, DataCostFunctor& costFunctor);
		
		// Get the unique edge description from key1 to key2 (essentially setting the reverse flag)
		// This function will fail if the edge is not unique
		EdgeDescription getUniqueEdgeDesc(const K& key1, const K& key2, extDirection parentDir);
		
		// return the number of edges a particular node has in the specified direction
		size_t getDegree(const K& key, extDirection dir);
		
		// return the number of vertices
		size_t getNumVertices() const { return m_vertexTable.size(); }
		
		// count the number of edges (SLOW)
		size_t countEdges() const;
		
		// print a node using ostream
		void printVertex(const K& key, bool printData = false) const;
		
		// validate the graph, looking for inconsistent links
		template<class Functor>
		void validate(Functor dataChecker);
		
		// Append nodes and merges nodes so that every vertex has 2 (or 0) edges in each direction
		template<class Functor>
		size_t removeTransitivity(Functor dataMerger);
		
		// attempt to merge two vertices along their shortest path with no guarentee they are linked
		template<class MergerFunctor>
		bool mergePath(const K& key1, const K& key2, extDirection parentDir, bool removeChild, bool usableChild, MergerFunctor dataMerger);
		
		// attempt to merge two vertices along their shortest path with no guarentee they are linked
		template<class DataCostFunctor, class MergerFunctor>
		bool mergeShortestPath(const K& key1, const K& key2, DataCostFunctor costFunctor, MergerFunctor dataMerger);
		
		// debug function to merge two vertices together
		template<class Functor>
		bool mergeWrapper(const K& key1, const K& key2, bool forceRemove, Functor dataMerger);
		
		// Iteratively visit each node
		template<class Functor>
		void iterativeVisit(Functor visitor);
		
		template<class Functor>
		void outputVertexConnectivity(Functor visitor) const;
		
	private:
	
		// Extract the shortest path between two vertices
		void extractShortestPath(VertexType* pSource, VertexType* pTarget, ShortestPathData& shortestPathData, KeyVec& path);
		
		// Run dijkstra's algorithm to find the shortest path between source and target using the cost functor specified
		template<class DataCostFunctor>
		void dijkstra(const K& sourceKey, ShortestPathData& shortestPathData, DataCostFunctor& costFunctor);
		
		//
		template<class DataCostFunctor>
		void greedyDirectedPath(const K& sourceKey, extDirection dir, KeySet& terminals, ShortestPathData& shortestPathData, DataCostFunctor& costFunctor);
		
		//
		template<class DataCostFunctor>		
		void ConstrainedDFS(VertexType* pCurrVertex, extDirection dir, VertexPtrSet vertexConstraints, 
										VertexPath currentPath, FeasiblePaths& solutions,
										size_t currLen, const size_t maxPathLen, DataCostFunctor& costFunctor);	
				
		//
		template<class DataCostFunctor>
		size_t getMinPathLength(const VertexPtrSet& vertexSet, DataCostFunctor costFunctor);											
		
		// Merge two vertices
		template<class Functor>
		bool merge(VertexType* parent, VertexType* child,  const extDirection parentsDir, const bool parentsReverse, bool removeChild, bool usableChild, Functor dataMerger);
		
		// Print a path of vertices
		template<class DataCostFunctor>
		void printVertexPath(const VertexPath& path, DataCostFunctor costFunctor);
		
		// calculate the length of a path
		template<class DataCostFunctor>
		size_t calculatePathLength(const VertexPath& path, DataCostFunctor costFunctor);
			
		VertexType* findVertex(const K& key) const;
		VertexTable m_vertexTable;
};

// the template implementation file
#include "DirectedGraphImpl.h"

#endif
