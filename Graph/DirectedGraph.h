#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H 1

#include "ContigNode.h"
#include "ContigPath.h"
#include <map>
#include <stdint.h>
#include <vector>

typedef std::vector<ContigPath> ContigPaths;

template<typename K, typename D>
struct Vertex
{
	typedef Vertex<K,D> VertexType;

	Vertex(const K& k, const D& d) : m_key(k), m_data(d) { }

	struct EdgeData
	{
		EdgeData(VertexType* node) : pVertex(node) { }

		VertexType* pVertex;

		bool operator==(const EdgeData& o) const
		{
			return pVertex == o.pVertex;
		}
	};

	/** Return the number of outgoing edges. */
	size_t numEdges() const
	{
		return m_edges.size();
	}

	/** Add an edge to this vertex. */
	void addEdge(VertexType* pNode);

	K m_key;
	D m_data;

	typedef typename std::vector<EdgeData> EdgeCollection;
	EdgeCollection m_edges;
};

typedef std::pair<ContigNode, unsigned> Constraint;
typedef std::vector<Constraint> Constraints;

template<typename D>
class DirectedGraph
{
	public:
		typedef ContigNode Node;
		typedef Vertex<Node, D> VertexType;

		/** Return the vertex specified by the given key. */
		const VertexType& operator[](const Node& key) const
		{
			const VertexType* pVertex = findVertex(key);
			assert(pVertex != NULL);
			return *pVertex;
		}

		/** Return the vertex specified by the given key. */
		VertexType& operator[](const Node& key)
		{
			VertexType* pVertex = findVertex(key);
			assert(pVertex != NULL);
			return *pVertex;
		}

		void addEdge(const Node& parent, const Node& child);
		void addVertex(const Node& key, const D& data);

		bool findSuperpaths(const Node& sourceKey,
				Constraints& constraints,
				ContigPaths& superPaths, unsigned& compCost) const;

		size_t getNumVertices() const { return m_vertexTable.size(); }
		size_t countEdges() const;
		void makeDistanceMap(const ContigPath& path,
				std::map<Node, int>& distanceMap) const;

	private:
		bool ConstrainedDFS(const VertexType* pCurrVertex,
				Constraints& constraints,
				Constraints::const_iterator nextConstraint,
				unsigned satisfied,
				ContigPath& path, ContigPaths& solutions,
				size_t currLen, unsigned& visitedCount) const;

		VertexType* findVertex(const Node& key)
		{
			return &m_vertexTable[key.index()];
		}

		const VertexType* findVertex(const Node& key) const
		{
			return &m_vertexTable[key.index()];
		}

		typedef typename std::vector<VertexType> VertexTable;
		VertexTable m_vertexTable;
};

#endif
