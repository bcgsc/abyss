#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H 1

#include "AffixIterator.h"
#include "ContigNode.h"
#include "ContigPath.h"
#include <iterator>
#include <map>
#include <ostream>
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

		friend std::ostream& operator <<(std::ostream& out,
				const EdgeData& o)
		{
			return out << o.node->m_key;
		}
	};

	/** Return the number of outgoing edges. */
	size_t out_degree() const
	{
		return m_edges.size();
	}

	/** Add an edge to this vertex. */
	void add_edge(VertexType* pNode);

	friend std::ostream& operator <<(std::ostream& out,
			const VertexType& o)
	{
		if (o.m_edges.empty())
			return out;
		out << '"' << o.m_key << "\" ->";
		if (o.m_edges.size() > 1)
			out << " {";
		std::copy(o.m_edges.begin(), o.m_edges.end(),
				affix_ostream_iterator<EdgeData>(out, " \"", "\""));
		if (o.m_edges.size() > 1)
			out << " }";
		return out;
	}

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

		/** Create an empty graph. */
		DirectedGraph() { }

		/** Create a graph with n vertices and zero edges. */
		DirectedGraph(unsigned n)
		{
			m_vertices.reserve(n);
			for (unsigned i = 0; i < n; ++i)
				add_vertex(Node(i));
			assert(m_vertices.size() == n);
		}

		/** Swap this graph with graph x. */
		void swap(DirectedGraph& x) { m_vertices.swap(x.m_vertices); }

		/** Return the vertex specified by the given key. */
		const VertexType& operator[](const Node& key) const
		{
			return m_vertices[key.index()];
		}

		/** Return the vertex specified by the given key. */
		VertexType& operator[](const Node& key)
		{
			return m_vertices[key.index()];
		}

		/** Remove all the edges and vertices from this graph. */
		void clear() { m_vertices.clear(); }

		void add_edge(const Node& parent, const Node& child);
		void add_vertex(const Node& key, const D& data = D());

		size_t num_vertices() const { return m_vertices.size(); }
		size_t num_edges() const;

		/** Returns the target vertex of edge e. */
		Node target(const typename VertexType::EdgeData& e) const
		{
			return e.node->m_key;
		}

		bool findSuperpaths(const Node& sourceKey,
				Constraints& constraints,
				ContigPaths& superPaths, unsigned& compCost) const;

		void makeDistanceMap(const ContigPath& path,
				std::map<Node, int>& distanceMap) const;

		friend std::ostream& operator <<(std::ostream& out,
				const DirectedGraph<D>& o)
		{
			std::copy(o.m_vertices.begin(), o.m_vertices.end(),
					std::ostream_iterator<VertexType>(out, "\n"));
			return out;
		}

	private:
		DirectedGraph(const DirectedGraph& x);
		DirectedGraph& operator =(const DirectedGraph& x);

		bool depthFirstSearch(const VertexType& currVertex,
				Constraints& constraints,
				Constraints::const_iterator nextConstraint,
				unsigned satisfied,
				ContigPath& path, ContigPaths& solutions,
				size_t currLen, unsigned& visitedCount) const;

		typedef typename std::vector<VertexType> VertexTable;
		VertexTable m_vertices;
};

#endif
