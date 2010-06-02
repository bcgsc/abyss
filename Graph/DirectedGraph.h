#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H 1

#include "ContigNode.h"
#include "ContigPath.h"
#include "Sense.h"
#include <map>
#include <stdint.h>
#include <vector>

typedef std::vector<ContigPath> ContigPaths;

typedef uint32_t LinearNumKey;

template<typename K, typename D>
struct Vertex
{
	typedef Vertex<K,D> VertexType;

	Vertex(const K& k, const D& d) : m_key(k), m_data(d) { }

	struct EdgeData
	{
		EdgeData(VertexType* node, bool rc)
			: pVertex(node), reverse(rc) { }

		VertexType* pVertex;
		bool reverse;

		bool operator==(const EdgeData& o) const
		{
			return pVertex == o.pVertex && reverse == o.reverse;
		}
	};

	/** Return the number of edges in the specified direction. */
	size_t numEdges(bool sense) const
	{
		return m_edges[sense].size();
	}

	// add an edge to the vertex in the specified direction
	void addEdge(VertexType* pNode, extDirection dir, bool reverse);

	K m_key;
	D m_data;

	typedef typename std::vector<EdgeData> EdgeCollection;
	EdgeCollection m_edges[NUM_DIRECTIONS];
};

typedef std::pair<ContigNode, unsigned> Constraint;
typedef std::vector<Constraint> Constraints;

template<typename D>
class DirectedGraph
{
	public:
		typedef Vertex<LinearNumKey, D> VertexType;

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

		void addEdge(const LinearNumKey& parent, extDirection dir,
				const ContigNode& child);
		void addVertex(const LinearNumKey& key, const D& data);

		bool findSuperpaths(const LinearNumKey& sourceKey,
				extDirection dir, Constraints& constraints,
				ContigPaths& superPaths, unsigned& compCost) const;

		size_t getNumVertices() const { return m_vertexTable.size(); }
		size_t countEdges() const;
		size_t calculatePathLength(const ContigPath& path) const;
		void makeDistanceMap(const ContigPath& path,
				std::map<ContigNode, int>& distanceMap) const;

	private:
		bool ConstrainedDFS(const VertexType* pCurrVertex,
				extDirection dir, Constraints& constraints,
				Constraints::const_iterator nextConstraint,
				unsigned satisfied,
				ContigPath& path, ContigPaths& solutions,
				size_t currLen, unsigned& visitedCount) const;

		VertexType* findVertex(const LinearNumKey& key);
		const VertexType* findVertex(const LinearNumKey& key) const;

		typedef typename std::vector<VertexType> VertexTable;
		VertexTable m_vertexTable;
};

#endif
