#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H 1

#include "AffixIterator.h"
#include "ContigNode.h"
#include <algorithm>
#include <cassert>
#include <iterator>
#include <ostream>
#include <utility>
#include <vector>

/** No properties. */
struct no_property
{
	friend std::ostream& operator <<(std::ostream& out,
			const no_property&)
	{
		return out;
	}

	friend std::istream& operator >>(std::istream& in, no_property&)
	{
		return in;
	}
};

/** A directed graph. */
template <typename VertexProp = no_property>
class DirectedGraph
{
  public:
	class Vertex;
	typedef typename std::vector<Vertex> Vertices;

/** Iterate through the vertices of this graph. */
class vertex_iterator {
	typedef typename Vertices::const_iterator const_iterator;

  public:
	vertex_iterator(const const_iterator& it) : m_it(it) { }
	const Vertex& operator *() const { return *m_it; }
	const Vertex* operator ->() const { return &*m_it; }

	bool operator !=(const vertex_iterator& it) const
	{
		return m_it != it.m_it;
	}

	vertex_iterator& operator ++() { ++m_it; return *this; }

  private:
	const_iterator m_it;
};

	class Edge;
	typedef typename std::vector<Edge> Edges;
	typedef typename Edges::const_iterator out_edge_iterator;

/** A vertex and its properties. */
class Vertex : public VertexProp
{
  public:
	Vertex() { }
	Vertex(const VertexProp& p) : VertexProp(p) { }

	/** Return an iterator to the edges of this vertex. */
	out_edge_iterator begin() const { return m_edges.begin(); }
	out_edge_iterator end() const { return m_edges.end(); }

	/** Return the first out edge of this vertex. */
	const Edge& front() const
	{
		assert(!m_edges.empty());
		return m_edges.front();
	}

	/** Return the number of outgoing edges. */
	size_t out_degree() const
	{
		return m_edges.size();
	}

	/** Add an edge to this vertex. */
	void add_edge(Vertex* v)
	{
		for (out_edge_iterator it = m_edges.begin();
				it != m_edges.end(); ++it)
			assert(v != &it->target());
		m_edges.push_back(v);
	}

	/** Remove the edge to v from this vertex. */
	void remove_edge(const Vertex& v)
	{
		m_edges.erase(find(m_edges.begin(), m_edges.end(), v));
	}

	/** Remove all out edges from this vertex. */
	void clear_out_edges()
	{
		m_edges.clear();
	}

	bool operator ==(const Vertex& v) const { return this == &v; }
	bool operator !=(const Vertex& v) const { return this != &v; }

  private:
	Edges m_edges;
};

/** A directed edge. */
class Edge
{
  public:
	Edge(Vertex* v) : m_target(v) { }

	/** Returns the target vertex of this edge. */
	const Vertex& target() const { return *m_target; }

	bool operator ==(const Vertex& v) { return m_target == &v; }

  private:
	/** The target vertex of this edge. */
	Vertex* m_target;
};

	public:
		typedef unsigned vertices_size_type;
		typedef unsigned degree_size_type;
		typedef ContigNode vertex_descriptor;
		typedef const Edge& edge_descriptor;

		/** Create an empty graph. */
		DirectedGraph() { }

		/** Create a graph with n vertices and zero edges. */
		DirectedGraph(vertices_size_type n) : m_vertices(n) { }

		/** Swap this graph with graph x. */
		void swap(DirectedGraph& x) { m_vertices.swap(x.m_vertices); }

		/** Return the vertex specified by the given descriptor. */
		const Vertex& operator[](vertex_descriptor v) const
		{
			return m_vertices[v.index()];
		}

		/** Return the vertex specified by the given descriptor. */
		Vertex& operator[](vertex_descriptor v)
		{
			return m_vertices[v.index()];
		}

		/** Returns an iterator-range to the vertices. */
		std::pair<vertex_iterator, vertex_iterator> vertices()
		{
			return make_pair(m_vertices.begin(), m_vertices.end());
		}

		/** Remove all the edges and vertices from this graph. */
		void clear() { m_vertices.clear(); }

		/** Add vertex v to the graph. */
		vertex_descriptor add_vertex(
				const VertexProp& data = VertexProp())
		{
			m_vertices.push_back(Vertex(data));
			return vertex(m_vertices.back());
		}

		/** Returns an iterator-range to the out edges of vertex u. */
		std::pair<out_edge_iterator, out_edge_iterator>
		out_edges(vertex_descriptor u)
		{
			const Vertex& v = (*this)[u];
			return make_pair(v.begin(), v.end());
		}

		/** Adds edge (u,v) to the graph. */
		void add_edge(vertex_descriptor u, vertex_descriptor v)
		{
			assert(u.index() < m_vertices.size());
			assert(v.index() < m_vertices.size());
			(*this)[u].add_edge(&(*this)[v]);
		}

		/** Remove the edge (u,v) from this graph. */
		void remove_edge(vertex_descriptor u, vertex_descriptor v)
		{
			(*this)[u].remove_edge((*this)[v]);
		}

		/** Remove all out edges from vertex v. */
		void clear_out_edges(vertex_descriptor v)
		{
			(*this)[v].clear_out_edges();
		}

		/** Remove vertex v from this graph. It is assumed that there
		 * are no edges to or from vertex v. It is best to call
		 * clear_vertex before remove_vertex.
		 */
		void remove_vertex(vertex_descriptor v)
		{
			unsigned i = v.index();
			if (i >= m_removed.size())
				m_removed.resize(i + 1);
			m_removed[i] = true;
		}

		/** Return the number of vertices. */
		size_t num_vertices() const { return m_vertices.size(); }

		/** Return an iterator to the vertex set of this graph. */
		vertex_iterator begin() const { return m_vertices.begin(); }
		vertex_iterator end() const { return m_vertices.end(); }

		/** Return the number of edges. */
		size_t num_edges() const
		{
			size_t n = 0;
			for (vertex_iterator it = m_vertices.begin();
					it != m_vertices.end(); ++it)
				n += it->out_degree();
			return n;
		}

		/** Return the out degree of the specified vertex. */
		degree_size_type out_degree(vertex_descriptor v) const
		{
			return (*this)[v].out_degree();
		}

		/** Return the nth vertex. */
		static vertex_descriptor vertex(vertices_size_type n)
		{
			return vertex_descriptor(n);
		}

		/** Return the descriptor of the specified vertex. */
		vertex_descriptor vertex(const Vertex& v) const
		{
			assert(&m_vertices[0] <= &v
					&& &v <= &m_vertices[0] + m_vertices.size());
			return vertex(&v - &m_vertices[0]);
		}

		/** Return the target vertex of the specified edge. */
		vertex_descriptor target(edge_descriptor e) const
		{
			return vertex(e.target());
		}

		friend std::ostream& operator <<(std::ostream& out,
				const DirectedGraph<VertexProp>& g)
		{
			for (vertex_iterator v = g.begin(); v != g.end(); ++v) {
				vertex_descriptor id = g.vertex(*v);
				if (g.is_removed(id))
					continue;
				if (sizeof (VertexProp) > 0)
					out << '"' << id << "\" [" << *v << "]\n";
				if (v->out_degree() == 0)
					continue;
				out << '"' << id << "\" ->";
				if (v->out_degree() > 1)
					out << " {";
				for (out_edge_iterator e = v->begin();
						e != v->end(); ++e)
					out << " \"" << g.target(*e) << '"';
				if (v->out_degree() > 1)
					out << " }";
				out << '\n';
			}
			return out;
		}

	protected:
		/** Return true if this vertex has been removed. */
		bool is_removed(vertex_descriptor v) const
		{
			unsigned i = v.index();
			return i < m_removed.size() ? m_removed[i] : false;
		}

	private:
		DirectedGraph(const DirectedGraph& x);
		DirectedGraph& operator =(const DirectedGraph& x);

		/** The set of vertices. */
		Vertices m_vertices;

		/** Flags indicating vertices that have been removed. */
		std::vector<bool> m_removed;
};

#endif
