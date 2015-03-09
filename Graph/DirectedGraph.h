#ifndef DIRECTEDGRAPH_H
#define DIRECTEDGRAPH_H 1

#include "Common/ContigNode.h"
#include "Graph/Properties.h"
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

/** A directed graph. */
template <typename VertexProp = no_property,
		 typename EdgeProp = no_property>
class DirectedGraph
{
	class Vertex;
	typedef typename std::vector<Vertex> Vertices;
	class Edge;
	typedef typename std::vector<Edge> Edges;

  public:
	// Graph
	typedef ContigNode vertex_descriptor;

	// IncidenceGraph
	typedef std::pair<vertex_descriptor, vertex_descriptor>
		edge_descriptor;
	typedef unsigned degree_size_type;

	// BidirectionalGraph
	typedef void in_edge_iterator;

	// VertexListGraph
	typedef unsigned vertices_size_type;

	// EdgeListGraph
	typedef unsigned edges_size_type;

	// PropertyGraph
	typedef VertexProp vertex_bundled;
	typedef VertexProp vertex_property_type;
	typedef EdgeProp edge_bundled;
	typedef EdgeProp edge_property_type;

	typedef boost::directed_tag directed_category;
	typedef boost::allow_parallel_edge_tag edge_parallel_category;
	struct traversal_category
		: boost::incidence_graph_tag,
		boost::adjacency_graph_tag,
		boost::vertex_list_graph_tag,
		boost::edge_list_graph_tag { };

/** Iterate through the vertices of this graph. */
class vertex_iterator
	: public std::iterator<std::input_iterator_tag,
		const vertex_descriptor>
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
	out_edge_iterator() { }
	out_edge_iterator(const const_iterator& it,
			vertex_descriptor src) : m_it(it), m_src(src) { }

	edge_descriptor operator *() const
	{
		return edge_descriptor(m_src, m_it->target());
	}

	bool operator ==(const out_edge_iterator& it) const
	{
		return m_it == it.m_it;
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

	const edge_property_type& get_property() const
	{
		return m_it->get_property();
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
	Vertex(const vertex_property_type& p) : m_prop(p) { }

	/** Return the properties of this vertex. */
	const vertex_property_type& get_property() const
	{
		return m_prop;
	}

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
	bool add_edge(vertex_descriptor v, const edge_property_type& ep)
	{
		m_edges.push_back(Edge(v, ep));
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

	/** Return the properties of the edge with target v. */
	edge_property_type& operator[](vertex_descriptor v)
	{
		typename Edges::iterator it
			= find(m_edges.begin(), m_edges.end(), v);
		assert(it != m_edges.end());
		return it->get_property();
	}

	/** Return the properties of the edge with target v. */
	const edge_property_type& operator[](vertex_descriptor v) const
	{
		typename Edges::const_iterator it
			= find(m_edges.begin(), m_edges.end(), v);
		assert(it != m_edges.end());
		return it->get_property();
	}

	/** Return true if edge (u,v) exists. */
	bool edge(vertex_descriptor v) const
	{
		return count(m_edges.begin(), m_edges.end(), v) > 0;
	}

	/** Remove edges that satisfy the predicate. */
	template <typename Predicate>
	void remove_edge_if(vertex_descriptor u, Predicate predicate)
	{
		typename Edges::iterator out = m_edges.begin();
		for (typename Edges::iterator it = m_edges.begin();
				it != m_edges.end(); ++it) {
			if (!predicate(edge_descriptor(u, it->target()))) {
				if (out != it)
					*out = *it;
				++out;
			}
		}
		m_edges.erase(out, m_edges.end());
	}

  private:
	Edges m_edges;
	vertex_property_type m_prop;
};

/** A directed edge. */
class Edge
{
  public:
	explicit Edge(vertex_descriptor v, const edge_property_type& ep)
		: m_target(v), m_ep(ep) { }

	/** Returns the target vertex of this edge. */
	vertex_descriptor target() const { return m_target; }

	/** Return true if the target of this edge is v. */
	bool operator ==(const vertex_descriptor& v) const
	{
		return m_target == v;
	}

	edge_property_type& get_property() { return m_ep; }
	const edge_property_type& get_property() const { return m_ep; }

  private:
	/** The target vertex of this edge. */
	vertex_descriptor m_target;
	edge_property_type m_ep;
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
	const vertex_property_type& operator[](vertex_descriptor u) const
	{
		vertices_size_type ui = get(vertex_index, *this, u);
		assert(ui < num_vertices());
		return m_vertices[ui].get_property();
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
	vertex_descriptor add_vertex(
			const vertex_property_type& vp = vertex_property_type())
	{
		m_vertices.push_back(Vertex(vp));
		return vertex_descriptor(num_vertices() - 1);
	}

	/** Returns an iterator-range to the out edges of vertex u. */
	std::pair<out_edge_iterator, out_edge_iterator>
	out_edges(vertex_descriptor u) const
	{
		vertices_size_type ui = get(vertex_index, *this, u);
		assert(ui < num_vertices());
		return m_vertices[ui].out_edges(u);
	}

	/** Returns an iterator-range to the adjacent vertices of
	 * vertex u. */
	std::pair<adjacency_iterator, adjacency_iterator>
	adjacent_vertices(vertex_descriptor u) const
	{
		vertices_size_type ui = get(vertex_index, *this, u);
		assert(ui < num_vertices());
		return m_vertices[ui].adjacent_vertices();
	}

	/** Adds edge (u,v) to this graph. */
	std::pair<edge_descriptor, bool>
	add_edge(vertex_descriptor u, vertex_descriptor v,
			const edge_property_type& ep = edge_property_type())
	{
		vertices_size_type ui = get(vertex_index, *this, u);
		assert(ui < num_vertices());
		assert(get(vertex_index, *this, v) < num_vertices());
		return make_pair(edge_descriptor(u, v),
				m_vertices[ui].add_edge(v, ep));
	}

	/** Remove the edge (u,v) from this graph. */
	void remove_edge(vertex_descriptor u, vertex_descriptor v)
	{
		vertices_size_type ui = get(vertex_index, *this, u);
		assert(ui < num_vertices());
		m_vertices[ui].remove_edge(v);
	}

	/** Remove the edge e from this graph. */
	void remove_edge(edge_descriptor e)
	{
		remove_edge(e.first, e.second);
	}

	/** Remove all out edges from vertex u. */
	void clear_out_edges(vertex_descriptor u)
	{
		vertices_size_type ui = get(vertex_index, *this, u);
		assert(ui < num_vertices());
		m_vertices[ui].clear_out_edges();
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

	/** Set the vertex_removed property. */
	void put(vertex_removed_t, vertex_descriptor u, bool flag)
	{
		vertices_size_type ui = get(vertex_index, *this, u);
		if (ui >= m_removed.size())
			m_removed.resize(ui + 1);
		m_removed[ui] = flag;
	}

	/** Remove vertex u from this graph. It is assumed that there
	 * are no edges to or from vertex u. It is best to call
	 * clear_vertex before remove_vertex.
	 */
	void remove_vertex(vertex_descriptor u)
	{
		put(vertex_removed, u, true);
	}

	/** Return the number of vertices. */
	vertices_size_type num_vertices() const
	{
		return m_vertices.size();
	}

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
		vertices_size_type ui = get(vertex_index, *this, u);
		assert(ui < num_vertices());
		return m_vertices[ui].out_degree();
	}

	/** Return the nth vertex. */
	static vertex_descriptor vertex(vertices_size_type n)
	{
		return vertex_descriptor(n);
	}

	/** Iterate through the edges of this graph. */
	std::pair<edge_iterator, edge_iterator> edges() const
	{
		std::pair<vertex_iterator, vertex_iterator> vit = vertices();
		return make_pair(edge_iterator(this, vit.first),
				edge_iterator(this, vit.second));
	}

	/** Return the edge (u,v) if it exists and a flag indicating
	 * whether the edge exists.
	 */
	std::pair<edge_descriptor, bool> edge(
			vertex_descriptor u, vertex_descriptor v) const
	{
		vertices_size_type ui = get(vertex_index, *this, u);
		assert(ui < num_vertices());
		return make_pair(edge_descriptor(u, v),
				m_vertices[ui].edge(v));
	}

	/** Return properties of edge e. */
	edge_property_type& operator[](edge_descriptor e)
	{
		vertices_size_type ui = get(vertex_index, *this, e.first);
		assert(ui < num_vertices());
		return m_vertices[ui][e.second];
	}

	/** Return properties of edge e. */
	const edge_property_type& operator[](edge_descriptor e) const
	{
		vertices_size_type ui = get(vertex_index, *this, e.first);
		assert(ui < num_vertices());
		return m_vertices[ui][e.second];
	}

	/** Remove edges that satisfy the predicate. */
	template <typename Predicate>
	void remove_edge_if(Predicate predicate)
	{
		unsigned i = 0;
		for (typename Vertices::iterator it = m_vertices.begin();
				it != m_vertices.end(); ++it)
			it->remove_edge_if(vertex(i++), predicate);
	}

	/** Return true if this vertex has been removed. */
	bool is_removed(vertex_descriptor u) const
	{
		vertices_size_type ui = get(vertex_index, *this, u);
		return ui < m_removed.size() ? m_removed[ui] : false;
	}

  private:
	DirectedGraph& operator =(const DirectedGraph& x);

	/** The set of vertices. */
	Vertices m_vertices;

	/** Flags indicating vertices that have been removed. */
	std::vector<bool> m_removed;
};

namespace std {
	template <typename VertexProp, typename EdgeProp>
	inline void swap(DirectedGraph<VertexProp, EdgeProp>& a,
			DirectedGraph<VertexProp, EdgeProp>& b) { a.swap(b); }
}

// IncidenceGraph

template <typename VP, typename EP>
std::pair<
	typename DirectedGraph<VP, EP>::out_edge_iterator,
	typename DirectedGraph<VP, EP>::out_edge_iterator>
out_edges(
		typename DirectedGraph<VP, EP>::vertex_descriptor u,
		const DirectedGraph<VP, EP>& g)
{
	return g.out_edges(u);
}

template <typename VP, typename EP>
typename DirectedGraph<VP, EP>::degree_size_type
out_degree(
		typename DirectedGraph<VP, EP>::vertex_descriptor u,
		const DirectedGraph<VP, EP>& g)
{
	return g.out_degree(u);
}

// AdjacencyGraph

template <typename VP, typename EP>
std::pair<
	typename DirectedGraph<VP, EP>::adjacency_iterator,
	typename DirectedGraph<VP, EP>::adjacency_iterator>
adjacent_vertices(
		typename DirectedGraph<VP, EP>::vertex_descriptor u,
		const DirectedGraph<VP, EP>& g)
{
	return g.adjacent_vertices(u);
}

// VertexListGraph

template <typename VP, typename EP>
typename DirectedGraph<VP, EP>::vertices_size_type
num_vertices(const DirectedGraph<VP, EP>& g)
{
	return g.num_vertices();
}

template <typename VP, typename EP>
typename DirectedGraph<VP, EP>::vertex_descriptor
vertex(typename DirectedGraph<VP, EP>::vertices_size_type ui, const DirectedGraph<VP, EP>& g)
{
	return g.vertex(ui);
}

template <typename VP, typename EP>
std::pair<typename DirectedGraph<VP, EP>::vertex_iterator,
	typename DirectedGraph<VP, EP>::vertex_iterator>
vertices(const DirectedGraph<VP, EP>& g)
{
	return g.vertices();
}

// EdgeListGraph

template <typename VP, typename EP>
typename DirectedGraph<VP, EP>::edges_size_type
num_edges(const DirectedGraph<VP, EP>& g)
{
	return g.num_edges();
}

template <typename VP, typename EP>
std::pair<typename DirectedGraph<VP, EP>::edge_iterator,
	typename DirectedGraph<VP, EP>::edge_iterator>
edges(const DirectedGraph<VP, EP>& g)
{
	return g.edges();
}

// AdjacencyMatrix

template <typename VP, typename EP>
std::pair<typename DirectedGraph<VP, EP>::edge_descriptor, bool>
edge(
	typename DirectedGraph<VP, EP>::vertex_descriptor u,
	typename DirectedGraph<VP, EP>::vertex_descriptor v,
	const DirectedGraph<VP, EP>& g)
{
	return g.edge(u, v);
}

// VertexMutableGraph

template <typename VP, typename EP>
typename DirectedGraph<VP, EP>::vertex_descriptor
add_vertex(DirectedGraph<VP, EP>& g)
{
	return g.add_vertex();
}

template <typename VP, typename EP>
void
remove_vertex(
		typename DirectedGraph<VP, EP>::vertex_descriptor u,
		DirectedGraph<VP, EP>& g)
{
	g.remove_vertex(u);
}

// EdgeMutableGraph

template <typename VP, typename EP>
void
clear_vertex(
		typename DirectedGraph<VP, EP>::vertex_descriptor u,
		DirectedGraph<VP, EP>& g)
{
	g.clear_vertex(u);
}

template <typename VP, typename EP>
std::pair<typename DirectedGraph<VP, EP>::edge_descriptor, bool>
add_edge(
	typename DirectedGraph<VP, EP>::vertex_descriptor u,
	typename DirectedGraph<VP, EP>::vertex_descriptor v,
	DirectedGraph<VP, EP>& g)
{
	return g.add_edge(u, v);
}

template <typename VP, typename EP>
void
remove_edge(
	typename DirectedGraph<VP, EP>::vertex_descriptor u,
	typename DirectedGraph<VP, EP>::vertex_descriptor v,
	DirectedGraph<VP, EP>& g)
{
	g.remove_edge(u, v);
}

template <typename VP, typename EP>
void
remove_edge(
		typename DirectedGraph<VP, EP>::edge_descriptor e,
		DirectedGraph<VP, EP>& g)
{
	g.remove_edge(e);
}

// MutableIncidenceGraph

template <typename VP, typename EP>
void
clear_out_edges(
		typename DirectedGraph<VP, EP>::vertex_descriptor u,
		DirectedGraph<VP, EP>& g)
{
	g.clear_out_edges(u);
}

// MutableEdgeListGraph

template <typename VP, typename EP, class Predicate>
void
remove_edge_if(Predicate predicate, DirectedGraph<VP, EP>& g)
{
	g.remove_edge_if(predicate);
}

// PropertyGraph

/** Return true if this vertex has been removed. */
template <typename VP, typename EP>
bool get(vertex_removed_t, const DirectedGraph<VP, EP>& g,
		typename DirectedGraph<VP, EP>::vertex_descriptor u)
{
	return g.is_removed(u);
}

template <typename VP, typename EP>
void put(vertex_removed_t tag, DirectedGraph<VP, EP>& g,
		typename DirectedGraph<VP, EP>::vertex_descriptor u,
		bool flag)
{
	g.put(tag, u, flag);
}

/** Return the properties of the edge of iterator eit. */
template <typename VP, typename EP>
const typename DirectedGraph<VP, EP>::edge_property_type&
get(edge_bundle_t, const DirectedGraph<VP, EP>&,
		typename DirectedGraph<VP, EP>::out_edge_iterator eit)
{
	return eit.get_property();
}

// PropertyGraph

template <typename VP, typename EP>
const VP&
get(vertex_bundle_t, const DirectedGraph<VP, EP>& g,
		typename DirectedGraph<VP, EP>::vertex_descriptor u)
{
	return g[u];
}

template <typename VP, typename EP>
const EP&
get(edge_bundle_t, const DirectedGraph<VP, EP>& g,
		typename DirectedGraph<VP, EP>::edge_descriptor e)
{
	return g[e];
}

// PropertyGraph vertex_index

namespace boost {
template <typename VP, typename EP>
struct property_map<DirectedGraph<VP, EP>, vertex_index_t>
{
	typedef ContigNodeIndexMap type;
	typedef type const_type;
};
}

template <typename VP, typename EP>
ContigNodeIndexMap
get(vertex_index_t, const DirectedGraph<VP, EP>&)
{
	return ContigNodeIndexMap();
}

template <typename VP, typename EP>
ContigNodeIndexMap::reference
get(vertex_index_t tag, const DirectedGraph<VP, EP>& g,
		typename DirectedGraph<VP, EP>::vertex_descriptor u)
{
	return get(get(tag, g), u);
}

// VertexMutablePropertyGraph

template <typename VP, typename EP>
typename DirectedGraph<VP, EP>::vertex_descriptor
add_vertex(const VP& vp, DirectedGraph<VP, EP>& g)
{
	return g.add_vertex(vp);
}

// EdgeMutablePropertyGraph

template <typename VP, typename EP>
std::pair<typename DirectedGraph<VP, EP>::edge_descriptor, bool>
add_edge(
	typename DirectedGraph<VP, EP>::vertex_descriptor u,
	typename DirectedGraph<VP, EP>::vertex_descriptor v,
	const typename DirectedGraph<VP, EP>::edge_property_type& ep,
	DirectedGraph<VP, EP>& g)
{
	return g.add_edge(u, v, ep);
}

// NamedGraph

template <typename VP, typename EP>
typename DirectedGraph<VP, EP>::vertex_descriptor
find_vertex(const std::string& name, const DirectedGraph<VP, EP>&)
{
	return find_vertex(name, g_contigNames);
}

template <typename VP, typename EP>
typename DirectedGraph<VP, EP>::vertex_descriptor
find_vertex(const std::string& name, bool sense,
		const DirectedGraph<VP, EP>&)
{
	return find_vertex(name, sense, g_contigNames);
}

#endif
