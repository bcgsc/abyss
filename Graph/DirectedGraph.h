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
#if HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP
# include <boost/graph/graph_traits.hpp>
#endif

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
	class Vertex;
	typedef typename std::vector<Vertex> Vertices;
	class Edge;
	typedef typename std::vector<Edge> Edges;
  public:
	typedef unsigned vertices_size_type;
	typedef ContigNode vertex_descriptor;
	typedef VertexProp vertex_property_type;
	typedef unsigned edges_size_type;
	typedef unsigned degree_size_type;
	typedef std::pair<vertex_descriptor, vertex_descriptor>
		edge_descriptor;
	typedef void in_edge_iterator;

#if HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP
	typedef boost::directed_tag directed_category;
	typedef boost::disallow_parallel_edge_tag edge_parallel_category;
	struct traversal_category
		: boost::incidence_graph_tag,
		boost::adjacency_graph_tag,
		boost::vertex_list_graph_tag,
		boost::edge_list_graph_tag { };
#else
	typedef void directed_category;
	typedef void edge_parallel_category;
	typedef void traversal_category;
#endif

/** Iterate through the vertices of this graph. */
class vertex_iterator
	: public std::iterator<std::input_iterator_tag, vertex_descriptor>
{
  public:
	vertex_iterator() { }
	explicit vertex_iterator(vertices_size_type v) : m_v(v) { }
	const vertex_descriptor& operator *() const { return m_v; }

	bool operator ==(const vertex_iterator& it) const
	{
		return m_v == it.m_v;
	}

	bool operator !=(const vertex_iterator& it) const
	{
		return m_v != it.m_v;
	}

	vertex_iterator& operator ++() { ++m_v; return *this; }
	vertex_iterator operator ++(int)
	{
		vertex_iterator it = *this;
		++*this;
		return it;
	}

  private:
	vertex_descriptor m_v;
};

/** Iterate through the out-edges. */
class out_edge_iterator
	: public std::iterator<std::input_iterator_tag, edge_descriptor>
{
	typedef typename Edges::const_iterator const_iterator;

  public:
	out_edge_iterator(const const_iterator& it,
			vertex_descriptor src) : m_it(it), m_src(src) { }

	edge_descriptor operator *() const
	{
		return edge_descriptor(m_src, m_it->target());
	}

	bool operator ==(const out_edge_iterator& it) const
	{
		return m_it != it.m_it;
	}

	bool operator !=(const out_edge_iterator& it) const
	{
		return m_it != it.m_it;
	}

	out_edge_iterator& operator ++() { ++m_it; return *this; }
	out_edge_iterator operator ++(int)
	{
		out_edge_iterator it = *this;
		++*this;
		return it;
	}

  private:
	const_iterator m_it;
	vertex_descriptor m_src;
};

/** Iterate through adjacent vertices. */
class adjacency_iterator : public Edges::const_iterator
{
	typedef typename Edges::const_iterator It;
  public:
	adjacency_iterator() { }
	adjacency_iterator(const It& it) : It(it) { }
	vertex_descriptor operator*() const
	{
		return It::operator*().target();
	}
};

/** Iterate through edges. */
class edge_iterator
	: public std::iterator<std::input_iterator_tag, edge_descriptor>
{
	void nextVertex()
	{
		vertex_iterator vlast = m_g->vertices().second;
		for (; m_vit != vlast; ++m_vit) {
			std::pair<adjacency_iterator, adjacency_iterator>
				adj = m_g->adjacent_vertices(*m_vit);
			if (adj.first != adj.second) {
				m_eit = adj.first;
				return;
			}
		}
		// Set m_eit to a known value.
		static const adjacency_iterator s_eitNULL;
		m_eit = s_eitNULL;
	}

  public:
	edge_iterator() { }
	edge_iterator(const DirectedGraph* g, const vertex_iterator& vit)
		: m_g(g), m_vit(vit)
	{
		nextVertex();
	}

	edge_descriptor operator*() const
	{
		return edge_descriptor(*m_vit, *m_eit);
	}

	bool operator==(const edge_iterator& it) const
	{
		return m_vit == it.m_vit && m_eit == it.m_eit;
	}

	bool operator!=(const edge_iterator& it) const
	{
		return !(*this == it);
	}

	edge_iterator& operator++()
	{
		if (++m_eit == m_g->adjacent_vertices(*m_vit).second) {
			++m_vit;
			nextVertex();
		}
		return *this;
	}

	edge_iterator operator++(int)
	{
		edge_iterator it = *this;
		++*this;
		return it;
	}

  private:
	const DirectedGraph* m_g;
	vertex_iterator m_vit;
	adjacency_iterator m_eit;
};

  private:
/** A vertex and its properties. */
class Vertex
{
  public:
	Vertex() { }
	Vertex(const VertexProp& p) : m_prop(p) { }

	/** Return the properties of this vertex. */
	const VertexProp& get_property() const { return m_prop; }

	/** Returns an iterator-range to the out edges of vertex u. */
	std::pair<out_edge_iterator, out_edge_iterator>
	out_edges(vertex_descriptor u) const
	{
		return make_pair(out_edge_iterator(m_edges.begin(), u),
				out_edge_iterator(m_edges.end(), u));
	}

	/** Returns an iterator-range to the adjacent vertices. */
	std::pair<adjacency_iterator, adjacency_iterator>
	adjacent_vertices() const
	{
		return make_pair(m_edges.begin(), m_edges.end());
	}

	/** Return the number of outgoing edges. */
	degree_size_type out_degree() const
	{
		return m_edges.size();
	}

	/** Add an edge to this vertex. */
	bool add_edge(vertex_descriptor v)
	{
		assert(count(m_edges.begin(), m_edges.end(), v) == 0);
		m_edges.push_back(Edge(v));
		return true;
	}

	/** Remove the edge to v from this vertex. */
	void remove_edge(vertex_descriptor v)
	{
		m_edges.erase(remove(m_edges.begin(), m_edges.end(), v),
				m_edges.end());
	}

	/** Remove all out edges from this vertex. */
	void clear_out_edges()
	{
		m_edges.clear();
	}

  private:
	Edges m_edges;
	VertexProp m_prop;
};

/** A directed edge. */
class Edge
{
  public:
	explicit Edge(vertex_descriptor v) : m_target(v) { }

	/** Returns the target vertex of this edge. */
	const vertex_descriptor& target() const { return m_target; }

	/** Return true if the target of this edge is v. */
	bool operator ==(const vertex_descriptor& v) const
	{
		return m_target == v;
	}

  private:
	/** The target vertex of this edge. */
	vertex_descriptor m_target;
};

  public:
	/** Create an empty graph. */
	DirectedGraph() { }

	/** Create a graph with n vertices and zero edges. */
	DirectedGraph(vertices_size_type n) : m_vertices(n) { }

	/** Swap this graph with graph x. */
	void swap(DirectedGraph& x)
	{
		m_vertices.swap(x.m_vertices);
		m_removed.swap(x.m_removed);
	}

	/** Return properties of vertex u. */
	const VertexProp& operator[](vertex_descriptor u) const
	{
		return m_vertices[u.index()].get_property();
	}

	/** Returns an iterator-range to the vertices. */
	std::pair<vertex_iterator, vertex_iterator> vertices() const
	{
		return make_pair(vertex_iterator(0),
			vertex_iterator(num_vertices()));
	}

	/** Remove all the edges and vertices from this graph. */
	void clear() { m_vertices.clear(); m_removed.clear(); }

	/** Add a vertex to this graph. */
	vertex_descriptor add_vertex(const VertexProp& vp = VertexProp())
	{
		m_vertices.push_back(Vertex(vp));
		return vertex_descriptor(num_vertices() - 1);
	}

	/** Returns an iterator-range to the out edges of vertex u. */
	std::pair<out_edge_iterator, out_edge_iterator>
	out_edges(vertex_descriptor u) const
	{
		assert(u.index() < num_vertices());
		return m_vertices[u.index()].out_edges(u);
	}

	/** Returns an iterator-range to the adjacent vertices of
	 * vertex u. */
	std::pair<adjacency_iterator, adjacency_iterator>
	adjacent_vertices(vertex_descriptor u) const
	{
		assert(u.index() < num_vertices());
		return m_vertices[u.index()].adjacent_vertices();
	}

	/** Adds edge (u,v) to this graph. */
	std::pair<edge_descriptor, bool>
	add_edge(vertex_descriptor u, vertex_descriptor v)
	{
		assert(u.index() < num_vertices());
		assert(v.index() < num_vertices());
		return make_pair(edge_descriptor(u, v),
				m_vertices[u.index()].add_edge(v));
	}

	/** Remove the edge (u,v) from this graph. */
	void remove_edge(vertex_descriptor u, vertex_descriptor v)
	{
		m_vertices[u.index()].remove_edge(v);
	}

	/** Remove the edge e from this graph. */
	void remove_edge(edge_descriptor e)
	{
		remove_edge(e.first, e.second);
	}

	/** Remove all out edges from vertex u. */
	void clear_out_edges(vertex_descriptor u)
	{
		m_vertices[u.index()].clear_out_edges();
	}

	/** Remove all edges to and from vertex u from this graph.
	 * O(V+E) */
	void clear_vertex(vertex_descriptor u)
	{
		clear_out_edges(u);
		std::pair<adjacency_iterator, adjacency_iterator>
			adj = adjacent_vertices(u);
		for (adjacency_iterator v = adj.first; v != adj.second; ++v)
			remove_edge(*v, u);
	}

	/** Remove vertex u from this graph. It is assumed that there
	 * are no edges to or from vertex u. It is best to call
	 * clear_vertex before remove_vertex.
	 */
	void remove_vertex(vertex_descriptor u)
	{
		unsigned i = u.index();
		if (i >= m_removed.size())
			m_removed.resize(i + 1);
		m_removed[i] = true;
	}

	/** Return the number of vertices. */
	vertices_size_type num_vertices() const { return m_vertices.size(); }

	/** Return the number of edges. */
	edges_size_type num_edges() const
	{
		edges_size_type n = 0;
		std::pair<vertex_iterator, vertex_iterator> vit = vertices();
		for (vertex_iterator v = vit.first; v != vit.second; ++v)
			n += out_degree(*v);
		return n;
	}

	/** Return the out degree of vertex u. */
	degree_size_type out_degree(vertex_descriptor u) const
	{
		return m_vertices[u.index()].out_degree();
	}

	/** Return the nth vertex. */
	static vertex_descriptor vertex(vertices_size_type n)
	{
		return vertex_descriptor(n);
	}

	/** Return the source vertex of the specified edge. */
	static vertex_descriptor source(edge_descriptor e)
	{
		return e.first;
	}

	/** Return the target vertex of the specified edge. */
	static vertex_descriptor target(edge_descriptor e)
	{
		return e.second;
	}

	/** Iterate through the edges of this graph. */
	std::pair<edge_iterator, edge_iterator> edges() const
	{
		std::pair<vertex_iterator, vertex_iterator> vit = vertices();
		return make_pair(edge_iterator(this, vit.first),
				edge_iterator(this, vit.second));
	}

	/** Return true if this vertex has been removed. */
	bool is_removed(vertex_descriptor u) const
	{
		unsigned i = u.index();
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

/** Output a GraphViz dot graph. */
template <typename Graph>
std::ostream& write_dot(std::ostream& out, const Graph& g)
{
	typedef typename Graph::vertex_iterator vertex_iterator;
	typedef typename Graph::adjacency_iterator adjacency_iterator;
	typedef typename Graph::vertex_property_type vertex_property_type;

	std::pair<vertex_iterator, vertex_iterator>
		vit = g.vertices();
	for (vertex_iterator u = vit.first; u != vit.second; ++u) {
		if (g.is_removed(*u))
			continue;
		if (sizeof (vertex_property_type) > 0)
			out << '"' << *u << "\" [" << g[*u] << "]\n";
		unsigned outdeg = g.out_degree(*u);
		if (outdeg == 0)
			continue;
		out << '"' << *u << "\" ->";
		if (outdeg > 1)
			out << " {";
		std::pair<adjacency_iterator, adjacency_iterator>
			adj = g.adjacent_vertices(*u);
		for (adjacency_iterator v = adj.first; v != adj.second; ++v)
			out << " \"" << *v << '"';
		if (outdeg > 1)
			out << " }";
		out << '\n';
	}
	return out;
}

/** Output a GraphViz dot graph. */
template <typename Graph>
struct dot_writer
{
	const Graph& g;
	dot_writer(const Graph& g) : g(g) { }
	friend std::ostream& operator<<(std::ostream& out,
			const dot_writer& o)
	{
		return write_dot<Graph>(out, o.g);
	}
};

template <typename VertexProp>
std::ostream& operator<<(std::ostream& out,
		const DirectedGraph<VertexProp>& g)
{
	return write_dot(out, g);
}

#endif
