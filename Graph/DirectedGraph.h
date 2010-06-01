#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H 1

#include "ContigNode.h"
#include "ContigPath.h"
#include "Sense.h"
#include <ostream>
#include <map>
#include <set>
#include <stdint.h>
#include <string>
#include <vector>

typedef std::vector<ContigPath> ContigPaths;

typedef uint32_t LinearNumKey;

enum VisitColor
{
	VC_WHITE,
	VC_BLACK
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

	/** Return the number of edges in the specified direction. */
	size_t numEdges(bool sense) const
	{
		return m_edges[sense].size();
	}

	// add an edge to the vertex in the specified direction
	void addEdge(VertexType* pNode, extDirection dir, bool reverse);
	
	// remove an edge
	void removeEdge(VertexType* pNode, extDirection dir, bool reverse);

	// check if the edge with the specified directionality exists
	bool edgeExists(const K& key, extDirection dir, bool reverse);

	// check if the described edge is unique
	bool isEdgeUnique(VertexType* pNode, extDirection dir, bool reverse);
	
	// get an edge
	EdgeCollectionIter getEdge(VertexType* pNode, extDirection dir, bool reverse, bool& found);
	
	// detect simple cycle
	bool detectSimpleCycle();

	K m_key;
	D m_data;	
	EdgeCollection m_edges[NUM_DIRECTIONS];
};

template<typename D>
class DirectedGraph 
{
	public:
		typedef typename std::vector<Vertex<LinearNumKey, D> > VertexTable;
		typedef typename VertexTable::iterator VertexTableIter;
		typedef typename VertexTable::const_iterator VertexTableConstIter;
		typedef Vertex<LinearNumKey, D> VertexType;

		typedef std::set<VertexType*> VertexCollection;
		typedef std::pair<LinearNumKey, VertexCollection> VertexComponent;
		typedef std::vector<VertexComponent> VertexComponentVector;

		typedef std::map<ContigNode, unsigned> KeyConstraintMap;
		typedef std::vector<LinearNumKey> KeyVec;

		typedef std::set<VertexType*> VertexPtrSet;

		typedef std::map<const VertexType*, size_t> DistanceMap;
		typedef std::map<const VertexType*, VisitColor> VisitedMap;
		typedef std::map<const VertexType*, const VertexType*> PreviousMap;
		typedef std::map<const VertexType*, extDirection> DirectionMap;

		struct ShortestPathData
		{
			DistanceMap distanceMap;
			VisitedMap visitedMap;
			PreviousMap previousMap;
			DirectionMap directionMap;
		};
	
		/** Return the vertex specified by the given key. */
		const VertexType& operator[](const LinearNumKey& key) const
		{
			const VertexType* pVertex = findVertex(key);
			assert(pVertex != NULL);
			return *pVertex;
		}

		/** Return the vertex specified by the given key. */
		VertexType& operator[](const LinearNumKey& key)
		{
			VertexType* pVertex = findVertex(key);
			assert(pVertex != NULL);
			return *pVertex;
		}

		/** Return the vertex data specified by the given key. */
		const D& getDataForVertex(const LinearNumKey& key) const
		{
			return (*this)[key].m_data;
		}

		void addEdge(const LinearNumKey& parent, extDirection dir,
				const ContigNode& child);
		void addVertex(const LinearNumKey& key, const D& data);
		void removeVertex(VertexType* pVertex);
		
		// reduce the graph with paired data		
		template<class ResolveFunctor>
		size_t reducePaired(ResolveFunctor& resolver);

		template<class ResolveFunctor, class DataMerger>
		bool attemptResolve(const LinearNumKey& key, extDirection dir,
				size_t maxCost, ResolveFunctor& resolver,
				DataMerger& merger);

		void generateComponents(VertexType* pVertex, extDirection dir,
				size_t maxCost, VertexComponentVector& outComponents);

		void accumulateVertices(VertexType* pVertex, extDirection dir,
				size_t currCost, size_t maxCost,
				VertexCollection& accumulator);

		bool findSuperpaths(const LinearNumKey& sourceKey,
				extDirection dir,
				const KeyConstraintMap& constraints,
				ContigPaths& superPaths, unsigned& compCost) const;

		// return the number of edges a particular node has in the specified direction
		size_t getDegree(const LinearNumKey& key, extDirection dir);
		
		// return the number of vertices
		size_t getNumVertices() const { return m_vertexTable.size(); }
		
		// count the number of edges (SLOW)
		size_t countEdges() const;
		
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
		template<class MergerFunctor>
		bool mergeShortestPath(const LinearNumKey& key1,
				const LinearNumKey& key2, MergerFunctor dataMerger);

		// debug function to merge two vertices together
		template<class Functor>
		bool mergeWrapper(const LinearNumKey& key1, const LinearNumKey& key2, bool forceRemove, Functor dataMerger);

		// Calculate the length of this path
		size_t calculatePathLength(const ContigPath& path) const;

		// Make a map of the distances to each node
		void makeDistanceMap(const ContigPath& path,
				std::map<ContigNode, int>& distanceMap) const;

	private:
		// Extract the shortest path between two vertices
		void extractShortestPath(const VertexType* pSource,
				const VertexType* pTarget,
				ShortestPathData& shortestPathData, KeyVec& path);

		// Run dijkstra's algorithm to find the shortest path between source and target using the cost functor specified
		void dijkstra(const LinearNumKey& sourceKey,
				ShortestPathData& shortestPathData);

		bool ConstrainedDFS(const VertexType* pCurrVertex,
				extDirection dir, const KeyConstraintMap& constraints,
				ContigPath& path, ContigPaths& solutions,
				size_t currLen, unsigned& visitedCount) const;

		size_t getMinPathLength(const VertexPtrSet& vertexSet);

		template<class Functor>
		bool merge(VertexType* parent, VertexType* child,  const extDirection parentsDir, const bool parentsReverse, bool removeChild, bool usableChild, Functor dataMerger);

		VertexType* findVertex(const LinearNumKey& key);
		const VertexType* findVertex(const LinearNumKey& key) const;
		VertexTable m_vertexTable;
};

#endif
