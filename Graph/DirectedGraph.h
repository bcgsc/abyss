#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H

#include "VisitAlgorithms.h"

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
	
	// find the direction that points to this key
	// the key HAS to be in the node
	size_t findEdgeIndex(const K& key);
	
	// add an edge to the vertex in the specified direction
	void addEdge(VertexType* pNode, extDirection dir, bool reverse);
	
	// remove an edge
	void removeEdge(VertexType* pNode, extDirection dir, bool reverse);
	
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
	typedef typename std::map<K, Vertex<K, D>* > VertexTable;
	typedef typename VertexTable::iterator VertexTableIter;
	typedef typename VertexTable::const_iterator VertexTableConstIter;
	typedef Vertex<K, D> VertexType;
	
	typedef std::set<VertexType*> VertexCollection;
	typedef std::pair<K, VertexCollection> VertexComponent;
	typedef std::vector<VertexComponent> VertexComponentVector;
	
	typedef std::vector<D*> DataCollection;
	typedef std::vector<DataCollection> DataComponents;
	
	public:
		DirectedGraph() { };
		~DirectedGraph();
		
		// add edge
		void addEdge(const K& parent, const K& child, extDirection dir, bool reverse);
		
		// add vertex
		VertexType* addVertex(const K& key, const D& data);
		
		// remove vertex
		template<class DataFunctor>
		void removeVertex(VertexType* pVertex, DataFunctor functor);
		
		// reduce the graph with paired data		
		template<class DataCostFunctor, class ResolveFunctor, class DataMerger>
		size_t reducePaired(size_t maxCost, DataCostFunctor& dataCost, ResolveFunctor& resolver, DataMerger& merger);
		
		// attempt to resolve 
		template<class DataCostFunctor, class ResolveFunctor, class DataMerger>
		bool attemptResolve(const K& key, extDirection dir, size_t maxCost, DataCostFunctor& dataCost, ResolveFunctor& resolver, DataMerger& merger);
		
		template<class DataCostFunctor>
		void generateComponents(VertexType* pVertex, extDirection dir, size_t maxCost, VertexComponentVector& outComponents, DataCostFunctor& dataCost);

		template<class DataCostFunctor>
		void accumulateVertices(VertexType* pVertex, extDirection dir, size_t currCost, size_t maxCost, VertexCollection& accumulator, DataCostFunctor& dataCost);
				
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
		
		// attempt to merge two vertices with no guarentee they are linked
		template<class Functor>
		bool mergeWrapper(const K& key1, const K& key2, Functor dataMerger);
		
		// Iteratively visit each node
		template<class Functor>
		void iterativeVisit(Functor visitor);
		
		template<class Functor>
		void bfsVisit(Functor f);
		
		void outputVertexConnectivity() const;
		
	private:
	
		// Merge two vertices
		template<class Functor>
		bool merge(VertexType* parent, VertexType* child,  const extDirection parentsDir, const bool parentsReverse, Functor dataMerger);
			
		VertexType* findVertex(const K& key) const;
		VertexTable m_vertexTable;
};

// the template implementation file
#include "DirectedGraphImpl.h"

#endif
