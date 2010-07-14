#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H 1

#include "AffixIterator.h"
#include "ContigNode.h"
#include <cassert>
#include <iterator>
#include <ostream>
#include <vector>

template<typename K, typename D>
class Vertex
{
  public:
	typedef Vertex<K,D> VertexType;
	class EdgeData;
	typedef typename std::vector<EdgeData> EdgeCollection;

	Vertex(const K& k, const D& d) : m_key(k), m_data(d) { }

	struct EdgeData
	{
		EdgeData(VertexType* node) : node(node) { }
		VertexType* node;

		/** Returns the target vertex of this edge. */
		K target() const
		{
			return node->m_key;
		}

		friend std::ostream& operator <<(std::ostream& out,
				const EdgeData& o)
		{
			return out << o.node->vertex();
		}
	};

	/** Return the key of this vertex. */
	const K& vertex() const { return m_key; }

	/** Return a collection of outgoing edges. */
	const EdgeCollection& out_edges() const { return m_edges; }

	/** Return the number of outgoing edges. */
	size_t out_degree() const
	{
		return m_edges.size();
	}

	/** Add an edge to this vertex. */
	void add_edge(VertexType* v)
	{
		for (typename EdgeCollection::const_iterator
				it = m_edges.begin(); it != m_edges.end(); ++it)
			assert(v != it->node);
		m_edges.push_back(v);
	}

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

  private:
	K m_key;
	D m_data;
	EdgeCollection m_edges;
};

template<typename D>
class DirectedGraph
{
	public:
		typedef ContigNode Node;
		typedef Vertex<Node, D> VertexType;
		typedef typename std::vector<VertexType> VertexTable;
		typedef typename VertexTable::const_iterator const_iterator;

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

		/** Adds vertex v to the graph. */
		void add_vertex(const Node& v, const D& data = D())
		{
			assert(m_vertices.size() == v.index());
			m_vertices.push_back(VertexType(v, data));
		}

		/** Adds edge (u,v) to the graph. */
		void add_edge(const Node& u, const Node& v)
		{
			assert(u.index() < m_vertices.size());
			assert(v.index() < m_vertices.size());
			(*this)[u].add_edge(&(*this)[v]);
		}

		/** Return the number of vertices. */
		size_t num_vertices() const { return m_vertices.size(); }

		/** Return an iterator to the vertex set of this graph. */
		const_iterator begin() const { return m_vertices.begin(); }
		const_iterator end() const { return m_vertices.end(); }

		/** Return the number of edges. */
		size_t num_edges() const
		{
			size_t n = 0;
			for (typename VertexTable::const_iterator it
					= m_vertices.begin(); it != m_vertices.end(); ++it)
				n += it->out_degree();
			return n;
		}

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

		VertexTable m_vertices;
};

#endif
