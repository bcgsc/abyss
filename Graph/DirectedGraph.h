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
		EdgeData(VertexType* node) : node(node) { }
		VertexType* node;
		bool operator==(const EdgeData& o) const
		{
			return node == o.node;
		}
	};

	/** Return the number of outgoing edges. */
	size_t out_degree() const
	{
		return m_edges.size();
	}

	/** Add an edge to this vertex. */
	void add_edge(VertexType* pNode);

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
			return m_vertexTable[key.index()];
		}

		/** Return the vertex specified by the given key. */
		VertexType& operator[](const Node& key)
		{
			return m_vertexTable[key.index()];
		}

		void add_edge(const Node& parent, const Node& child);
		void add_vertex(const Node& key, const D& data = D());

		size_t num_vertices() const { return m_vertexTable.size(); }
		size_t num_edges() const;

		bool findSuperpaths(const Node& sourceKey,
				Constraints& constraints,
				ContigPaths& superPaths, unsigned& compCost) const;

		void makeDistanceMap(const ContigPath& path,
				std::map<Node, int>& distanceMap) const;

	private:
		bool depthFirstSearch(const VertexType& currVertex,
				Constraints& constraints,
				Constraints::const_iterator nextConstraint,
				unsigned satisfied,
				ContigPath& path, ContigPaths& solutions,
				size_t currLen, unsigned& visitedCount) const;

		typedef typename std::vector<VertexType> VertexTable;
		VertexTable m_vertexTable;
};

#endif
