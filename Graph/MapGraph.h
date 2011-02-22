#ifndef MAPGRAPH_H
#define MAPGRAPH_H 1

#include "Graph.h"
#include <cassert>
#include <map>
#include <utility> // for make_pair

/**
 * A std::map is a model of a directed graph.
 * @author Shaun Jackman <sjackman@bcgsc.ca>
 */
template <typename V, typename T>
struct graph_traits<std::map<V, T> > {
	// Graph
	typedef V vertex_descriptor;
#if HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP
	typedef boost::directed_tag directed_category;
	struct traversal_category
		: boost::adjacency_graph_tag,
		boost::vertex_list_graph_tag
		{ };
	typedef boost::disallow_parallel_edge_tag edge_parallel_category;
#else
	typedef void directed_category;
	typedef void traversal_category;
	typedef void edge_parallel_category;
#endif

	// IncidenceGraph
	typedef std::pair<vertex_descriptor, vertex_descriptor>
		edge_descriptor;
	typedef unsigned degree_size_type;

	/** Iterate through the out edges of a vertex. */
	struct out_edge_iterator
		: std::map<V, typename T::mapped_type>::const_iterator
	{
		typedef typename std::map<V,
				typename T::mapped_type>::const_iterator It;
		out_edge_iterator(const It& it, vertex_descriptor src)
			: It(it), m_src(src) { }
		edge_descriptor operator*() const
		{
			return edge_descriptor(m_src, It::operator*().first);
		}
	  private:
		vertex_descriptor m_src;
	};

	// BidirectionalGraph
	typedef void in_edge_iterator;

	// AdjacencyGraph

	/** Iterate through the adjacent vertices of a vertex. */
	struct adjacency_iterator
		: std::map<V, typename T::mapped_type>::const_iterator
	{
		typedef typename std::map<V,
				typename T::mapped_type>::const_iterator It;
		adjacency_iterator(const It& it) : It(it) { }
		const vertex_descriptor& operator*() const
		{
			return It::operator*().first;
		}
	};

	// VertexListGraph

	/** Iterate through the vertices of this graph. */
	struct vertex_iterator : std::map<V, T>::const_iterator
	{
		typedef typename std::map<V, T>::const_iterator It;
		vertex_iterator(const It& it) : It(it) { }
		const vertex_descriptor& operator*() const
		{
			return It::operator*().first;
		}
	};

	typedef unsigned vertices_size_type;

	// EdgeListGraph
	typedef void edge_iterator;
	typedef void edges_size_type;
};

// IncidenceGraph

template <typename V, typename T>
std::pair<
	typename graph_traits<std::map<V, T> >::out_edge_iterator,
	typename graph_traits<std::map<V, T> >::out_edge_iterator>
out_edges(
		typename graph_traits<std::map<V, T> >::vertex_descriptor u,
		const std::map<V, T>& g)
{
	typedef typename graph_traits<std::map<V, T> >::out_edge_iterator
		out_edge_iterator;
	typename std::map<V, T>::const_iterator it = g.find(u);
	assert(it != g.end());
	return make_pair(
		out_edge_iterator(it->second.begin(), u),
		out_edge_iterator(it->second.end(), u));
}

template <typename V, typename T>
typename graph_traits<std::map<V, T> >::degree_size_type
out_degree(
		typename graph_traits<std::map<V, T> >::vertex_descriptor u,
		const std::map<V, T>& g)
{
	typename std::map<V, T>::const_iterator it = g.find(u);
	assert(it != g.end());
	return it->second.size();
}

// BidirectionalGraph

template <typename V, typename T>
typename graph_traits<std::map<V, T> >::degree_size_type
in_degree(
		typename graph_traits<std::map<V, T> >::vertex_descriptor u,
		const std::map<V, T>& g)
{
	return out_degree(~u, g);
}

// AdjacencyGraph

template <typename V, typename T>
std::pair<typename graph_traits<std::map<V, T> >::adjacency_iterator,
	typename graph_traits<std::map<V, T> >::adjacency_iterator>
adjacent_vertices(
		typename graph_traits<std::map<V, T> >::vertex_descriptor u,
		const std::map<V, T>& g)
{
	typename std::map<V, T>::const_iterator it = g.find(u);
	assert(it != g.end());
	return make_pair(it->second.begin(), it->second.end());
}

// VertexListGraph

template <typename V, typename T>
std::pair<typename graph_traits<std::map<V, T> >::vertex_iterator,
	typename graph_traits<std::map<V, T> >::vertex_iterator>
vertices(const std::map<V, T>& g)
{
	return make_pair(g.begin(), g.end());
}

// VertexMutableGraph

template <typename V, typename T>
void remove_vertex(
		typename graph_traits<std::map<V, T> >::vertex_descriptor u,
		std::map<V, T>& g)
{
	unsigned n = g.erase(u);
	assert(n == 1);
	(void)n;
}

// EdgeMutableGraph

template <typename V, typename T>
void clear_out_edges(
		typename graph_traits<std::map<V, T> >::vertex_descriptor u,
		std::map<V, T>& g)
{
	typename std::map<V, T>::iterator it = g.find(u);
	assert(it != g.end());
	it->second.clear();
}

template <typename V, typename T>
std::pair<
	typename graph_traits<std::map<V, T> >::edge_descriptor, bool>
add_edge(
		typename graph_traits<std::map<V, T> >::vertex_descriptor u,
		typename graph_traits<std::map<V, T> >::vertex_descriptor v,
		std::map<V, T>& g)
{
	typename edge_property<std::map<V, T> >::type ep;
	return make_pair(std::make_pair(u, v),
			g[u].insert(make_pair(v, ep)).second);
}

template <typename V, typename T>
void remove_edge(
		typename graph_traits<std::map<V, T> >::vertex_descriptor u,
		typename graph_traits<std::map<V, T> >::vertex_descriptor v,
		std::map<V, T>& g)
{
	unsigned n = g[u].erase(v);
	assert(n == 1);
	(void)n;
}

// PropertyGraph

template <typename V, typename T>
no_property get(vertex_bundle_t, const std::map<V, T>&,
		typename graph_traits<std::map<V, T> >::vertex_descriptor)
{
	return no_property();
}

template <typename V, typename T>
bool get(vertex_removed_t, const std::map<V, T>&,
		typename graph_traits<std::map<V, T> >::vertex_descriptor)
{
	return false;
}

template <typename V, typename T>
const typename T::mapped_type&
get(edge_bundle_t, const std::map<V, T>& g,
		typename graph_traits<std::map<V, T> >::edge_descriptor e)
{
	typename std::map<V, T>::const_iterator u = g.find(source(e, g));
	assert(u != g.end());
	typename T::const_iterator v = u->second.find(target(e, g));
	assert(v != u->second.end());
	return v->second;
}

template <typename V, typename T>
const typename T::mapped_type&
get(edge_bundle_t, const std::map<V, T>&,
		typename graph_traits<std::map<V, T> >::out_edge_iterator eit)
{
	return eit->second;
}

// VertexMutablePropertyGraph

template <typename V, typename T>
class vertex_property<std::map<V, T> > {
  public:
	typedef no_property type;
};

// EdgeMutablePropertyGraph

template <typename V, typename T>
class edge_property<std::map<V, T> > {
  public:
	typedef typename T::mapped_type type;
};

template <typename V, typename T>
std::pair<
	typename graph_traits<std::map<V, T> >::edge_descriptor, bool>
add_edge(
		typename graph_traits<std::map<V, T> >::vertex_descriptor u,
		typename graph_traits<std::map<V, T> >::vertex_descriptor v,
		typename edge_property<std::map<V, T> >::type ep,
		std::map<V, T>& g)
{
	return make_pair(std::make_pair(u, v),
			g[u].insert(std::make_pair(v, ep)).second);
}

#endif
