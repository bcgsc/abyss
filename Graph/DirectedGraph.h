#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H 1

#include "AffixIterator.h"
#include "ContigNode.h"
#include <cassert>
#include <iterator>
#include <ostream>
#include <vector>

template<typename D>
class Vertex
{
  public:
	typedef Vertex<D> VertexType;
	class Edge;
	typedef typename std::vector<Edge> Edges;

	Vertex() { }
	Vertex(const D& d) : m_data(d) { }

	class Edge
	{
	  public:
		Edge(VertexType* v) : m_target(v) { }

		/** Returns the target vertex of this edge. */
		const VertexType& target() const { return *m_target; }

		friend std::ostream& operator <<(std::ostream& out,
				const Edge& e)
		{
			return out << e.m_target;
		}

	  private:
		/** The target vertex of this edge. */
		VertexType* m_target;
	};

	/** Return a collection of outgoing edges. */
	const Edges& out_edges() const { return m_edges; }

	/** Return the number of outgoing edges. */
	size_t out_degree() const
	{
		return m_edges.size();
	}

	/** Add an edge to this vertex. */
	void add_edge(VertexType* v)
	{
		for (typename Edges::const_iterator it = m_edges.begin();
				it != m_edges.end(); ++it)
			assert(v != &it->target());
		m_edges.push_back(v);
	}

	bool operator ==(const VertexType& v) const { return this == &v; }
	bool operator !=(const VertexType& v) const { return this != &v; }

	friend std::ostream& operator <<(std::ostream& out,
			const VertexType& o)
	{
		if (o.m_edges.empty())
			return out;
		out << '"' << &o << "\" ->";
		if (o.m_edges.size() > 1)
			out << " {";
		std::copy(o.m_edges.begin(), o.m_edges.end(),
				affix_ostream_iterator<Edge>(out, " \"", "\""));
		if (o.m_edges.size() > 1)
			out << " }";
		return out;
	}

  private:
	D m_data;
	Edges m_edges;
};

template<typename D>
class DirectedGraph
{
	public:
		typedef ContigNode Node;
		typedef Vertex<D> VertexType;
		typedef typename std::vector<VertexType> Vertices;
		typedef typename Vertices::const_iterator const_iterator;
		typedef typename VertexType::Edge Edge;
		typedef typename VertexType::Edges Edges;

		/** Create an empty graph. */
		DirectedGraph() { }

		/** Create a graph with n vertices and zero edges. */
		DirectedGraph(unsigned n) : m_vertices(n) { }

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
			m_vertices.push_back(VertexType(data));
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
			for (typename Vertices::const_iterator it
					= m_vertices.begin(); it != m_vertices.end(); ++it)
				n += it->out_degree();
			return n;
		}

		/** Return the out degree of the specified vertex. */
		unsigned out_degree(const Node& v) const
		{
			return (*this)[v].out_degree();
		}

		/** Return the in degree of the specified vertex. */
		unsigned in_degree(const Node& v) const
		{
			return (*this)[~v].out_degree();
		}

		/** Return the in degree of the specified vertex. */
		unsigned in_degree(const VertexType& v) const
		{
			return in_degree(vertex(v));
		}

		/** Return the nth vertex. */
		static Node vertex(unsigned n)
		{
			return Node(n);
		}

		/** Return the descriptor of the specified vertex. */
		Node vertex(const VertexType& v) const
		{
			assert(&m_vertices[0] <= &v
					&& &v <= &m_vertices[0] + m_vertices.size());
			return vertex(&v - &m_vertices[0]);
		}

		/** Return the target vertex of the specified edge. */
		Node target(const Edge& e) const
		{
			return vertex(e.target());
		}

		friend std::ostream& operator <<(std::ostream& out,
				const DirectedGraph<D>& g)
		{
			for (const_iterator v = g.begin(); v != g.end(); ++v) {
				if (v->out_degree() == 0)
					continue;
				out << '"' << g.vertex(*v) << "\" ->";
				if (v->out_degree() > 1)
					out << " {";
				for (typename Edges::const_iterator e
						= v->out_edges().begin();
						e != v->out_edges().end(); ++e)
					out << " \"" << g.target(*e) << '"';
				if (v->out_degree() > 1)
					out << " }";
				out << '\n';
			}
			return out;
		}

	private:
		DirectedGraph(const DirectedGraph& x);
		DirectedGraph& operator =(const DirectedGraph& x);

		Vertices m_vertices;
};

#endif
