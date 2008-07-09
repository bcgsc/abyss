#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H

#include "Timer.h"
#include "SetOperations.h"
#include <queue>
#include <list>
#include <google/sparse_hash_map>
#include <ext/hash_map>
#include <iostream>
#include <string>


typedef uint32_t LinearNumKey;

enum VisitColor
{
	VC_WHITE,
	VC_GRAY,
	VC_BLACK
};

// Work around gnu_cxx's lack of support for hashing strings
namespace __gnu_cxx                                                                                 
{                                                                                             
  template<> struct hash< std::string >
  {
    size_t operator()( const std::string& x ) const                                           
    {                                                                                         
      return hash< const char* >()( x.c_str() );
    }                                                                                         
  };                                                                                          
}          

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
		
	// typedefs for convienance
	typedef Vertex<K,D> VertexType;
	
	struct EdgeData
	{
		VertexType* pVertex;
		bool reverse;
		
		bool operator==(const EdgeData& e2) const
		{
			return (this->pVertex == e2.pVertex && this->reverse == e2.reverse);
		}		
	};
	

		
	typedef typename std::vector<EdgeData> EdgeCollection;
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
		
	K m_key;
	D m_data;	
	EdgeCollection m_edges[NUM_DIRECTIONS];
	
};

template<typename D>
class DirectedGraph 
{
	public:
		// Lots of typedefs
		typedef typename std::vector<Vertex<LinearNumKey, D>* > VertexTable;
		//typedef typename std::map<K, Vertex<K, D>* > VertexTable;
		typedef typename VertexTable::iterator VertexTableIter;
		typedef typename VertexTable::const_iterator VertexTableConstIter;
		typedef Vertex<LinearNumKey, D> VertexType;
		
		typedef std::set<VertexType*> VertexCollection;
		typedef std::pair<LinearNumKey, VertexCollection> VertexComponent;
		typedef std::vector<VertexComponent> VertexComponentVector;
		
		typedef std::set<LinearNumKey> KeySet;
		typedef std::map<LinearNumKey, int> KeyIntMap;
		typedef std::vector<LinearNumKey> KeyVec;

		typedef std::set<VertexType*> VertexPtrSet;
		
		struct PathNode
		{
			LinearNumKey key;
			bool isRC;
			
			friend std::ostream& operator<<(std::ostream& out, const PathNode& object)
			{
				out << "(" << object.key << "," << object.isRC << ")";
				return out;
			}
		};
		
		typedef std::vector<PathNode> VertexPath;
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
	
		DirectedGraph() { };
		DirectedGraph(const size_t sizeHint) { m_vertexTable.reserve(sizeHint); }
		~DirectedGraph();
		
		// get the data pointer
		const D& getDataForVertex(const LinearNumKey& key) const { VertexType* pVertex = findVertex(key); assert(pVertex != NULL); return pVertex->m_data; }
		
		// add edge
		void addEdge(const LinearNumKey& parent, const LinearNumKey& child, extDirection dir, bool reverse);
		
		// add vertex
		VertexType* addVertex(const LinearNumKey& key, const D& data);
		
		// remove vertex
		void removeVertex(VertexType* pVertex);
		
		// reduce the graph with paired data		
		template<class ResolveFunctor>
		size_t reducePaired(ResolveFunctor& resolver);
		
		// attempt to resolve 
		template<class DataCostFunctor, class ResolveFunctor, class DataMerger>
		bool attemptResolve(const LinearNumKey& key, extDirection dir, size_t maxCost, DataCostFunctor& dataCost, ResolveFunctor& resolver, DataMerger& merger);
		
		template<class DataCostFunctor>
		void generateComponents(VertexType* pVertex, extDirection dir, size_t maxCost, VertexComponentVector& outComponents, DataCostFunctor& dataCost);

		template<class DataCostFunctor>
		void accumulateVertices(VertexType* pVertex, extDirection dir, size_t currCost, size_t maxCost, VertexCollection& accumulator, DataCostFunctor& dataCost);
		
		template<class DataCostFunctor>
		bool findSuperpaths(const LinearNumKey& sourceKey, extDirection dir, const KeyIntMap& keyConstraints, FeasiblePaths& superPaths, 
				DataCostFunctor& costFunctor, int maxNumPaths, int maxCompCost, int& compCost);
		
		// Get the unique edge description from key1 to key2 (essentially setting the reverse flag)
		// This function will fail if the edge is not unique
		EdgeDescription getUniqueEdgeDesc(const LinearNumKey& key1, const LinearNumKey& key2, extDirection parentDir);
		
		// return the number of edges a particular node has in the specified direction
		size_t getDegree(const LinearNumKey& key, extDirection dir);
		
		// return the number of vertices
		size_t getNumVertices() const { return m_vertexTable.size(); }
		
		// count the number of edges (SLOW)
		size_t countEdges() const;
		
		// print a node using ostream
		void printVertex(const LinearNumKey& key) const;
		
		// validate the graph, looking for inconsistent links
		template<class Functor>
		void validate(Functor dataChecker);
		
		// Append nodes and merges nodes so that every vertex has 2 (or 0) edges in each direction
		template<class Functor>
		size_t removeTransitivity(Functor dataMerger);
		
		// attempt to merge two vertices along their shortest path with no guarentee they are linked
		template<class MergerFunctor>
		bool mergePath(const LinearNumKey& key1, const LinearNumKey& key2, extDirection parentDir, bool removeChild, bool usableChild, MergerFunctor dataMerger);
		
		// attempt to merge two vertices along their shortest path with no guarentee they are linked
		template<class DataCostFunctor, class MergerFunctor>
		bool mergeShortestPath(const LinearNumKey& key1, const LinearNumKey& key2, DataCostFunctor costFunctor, MergerFunctor dataMerger);
		
		// debug function to merge two vertices together
		template<class Functor>
		bool mergeWrapper(const LinearNumKey& key1, const LinearNumKey& key2, bool forceRemove, Functor dataMerger);
		
		// Calculate the length of this path
		template<typename DataCostFunctor>
		size_t calculatePathLength(const VertexPath& path, DataCostFunctor costFunctor);
		
		// Make a map of the distances to each node
		template<typename DataCostFunctor>
		void makeDistanceMap(const VertexPath& path, DataCostFunctor costFunctor, std::map<LinearNumKey, int>& distanceMap);
		
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
		void dijkstra(const LinearNumKey& sourceKey, ShortestPathData& shortestPathData, DataCostFunctor& costFunctor);
		
		//
		template<class DataCostFunctor>
		void greedyDirectedPath(const LinearNumKey& sourceKey, extDirection dir, KeySet& terminals, ShortestPathData& shortestPathData, DataCostFunctor& costFunctor);
		
		//
		template<class DataCostFunctor>		
		void ConstrainedDFS(VertexType* pCurrVertex, extDirection dir, const KeyIntMap keyConstraints, 
										VertexPath currentPath, FeasiblePaths& solutions,
										size_t currLen, DataCostFunctor& costFunctor, int maxNumPaths, int maxCompCost);	
				
		//
		template<class DataCostFunctor>
		size_t getMinPathLength(const VertexPtrSet& vertexSet, DataCostFunctor costFunctor);											
		
		// Merge two vertices
		template<class Functor>
		bool merge(VertexType* parent, VertexType* child,  const extDirection parentsDir, const bool parentsReverse, bool removeChild, bool usableChild, Functor dataMerger);
		
			
		VertexType* findVertex(const LinearNumKey& key) const;
		VertexTable m_vertexTable;
};

// the template implementation file
#include "DirectedGraphImpl.h"

#endif
