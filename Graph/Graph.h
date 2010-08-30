#ifndef GRAPH_H
#define GRAPH_H 1

#include <utility> // for pair

// Graph

template <typename G>
struct graph_traits {
	// Graph
	typedef typename G::vertex_descriptor vertex_descriptor;
	typedef typename G::directed_category directed_category;
	typedef typename G::traversal_category traversal_category;
	typedef typename G::edge_parallel_category edge_parallel_category;

	// IncidenceGraph
	typedef typename G::edge_descriptor edge_descriptor;
	typedef typename G::out_edge_iterator out_edge_iterator;
	typedef typename G::degree_size_type degree_size_type;

	// BidirectionalGraph
	typedef typename G::in_edge_iterator in_edge_iterator;

	// AdjacencyGraph
	typedef typename G::adjacency_iterator adjacency_iterator;

	// VertexListGraph
	typedef typename G::vertex_iterator vertex_iterator;
	typedef typename G::vertices_size_type vertices_size_type;

	// EdgeListGraph
	typedef typename G::edge_iterator edge_iterator;
	typedef typename G::edges_size_type edges_size_type;
};

// IncidenceGraph

template <class G>
std::pair<
	typename G::out_edge_iterator,
	typename G::out_edge_iterator>
out_edges(typename G::vertex_descriptor u, const G& g)
{
	return g.out_edges(u);
}

#if !HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP
template <class G>
typename graph_traits<G>::vertex_descriptor
source(std::pair<typename graph_traits<G>::vertex_descriptor,
		typename graph_traits<G>::vertex_descriptor> e, const G&)
{
	return e.first;
}

template <class G>
typename graph_traits<G>::vertex_descriptor
target(std::pair<typename graph_traits<G>::vertex_descriptor,
		typename graph_traits<G>::vertex_descriptor> e, const G&)
{
	return e.second;
}
#endif

template <class G>
typename G::degree_size_type
out_degree(typename G::vertex_descriptor u, const G& g)
{
	return g.out_degree(u);
}

// BidirectionalGraph

template <class G>
typename G::degree_size_type
in_degree(typename G::vertex_descriptor u, const G& g)
{
	return g.in_degree(u);
}

// AdjacencyGraph

template <class G>
std::pair<
	typename G::adjacency_iterator,
	typename G::adjacency_iterator>
adjacent_vertices(typename G::vertex_descriptor u, const G& g)
{
	return g.adjacent_vertices(u);
}

// VertexListGraph

template <class G>
typename G::vertices_size_type
num_vertices(const G& g)
{
	return g.num_vertices();
}

template <class G>
std::pair<typename G::vertex_iterator, typename G::vertex_iterator>
vertices(const G& g)
{
	return g.vertices();
}

// EdgeListGraph

template <class G>
typename G::edges_size_type
num_edges(const G& g)
{
	return g.num_edges();
}

template <class G>
std::pair<typename G::edge_iterator, typename G::edge_iterator>
edges(const G& g)
{
	return g.edges();
}

// VertexMutableGraph

template <class G>
typename G::vertex_descriptor
add_vertex(G& g)
{
	return g.add_vertex();
}

template <class G>
void
clear_vertex(typename G::vertex_descriptor u, G& g)
{
	g.clear_vertex(u);
}

template <class G>
void
remove_vertex(typename G::vertex_descriptor u, G& g)
{
	g.remove_vertex(u);
}

// EdgeMutableGraph

template <class G>
std::pair<typename G::edge_descriptor, bool>
add_edge(
	typename G::vertex_descriptor u,
	typename G::vertex_descriptor v,
	G& g)
{
	return g.add_edge(u, v);
}

template <class G>
void
remove_edge(
	typename G::vertex_descriptor u,
	typename G::vertex_descriptor v,
	G& g)
{
	return g.remove_edge(u, v);
}

template <class G>
void
remove_edge(typename G::edge_descriptor e, G& g)
{
	g.remove_edge(e);
}

// Properties

/** A vertex bundle property. */
enum vertex_bundle_t { vertex_bundle };

/** A property indicating that this vertex has been removed. */
enum vertex_removed_t { vertex_removed };

/** An edge bundle property. */
enum edge_bundle_t { edge_bundle };

/** The distance between two vertices. */
enum edge_distance_t { edge_distance };

// VertexMutablePropertyGraph

template <class Graph>
struct vertex_property {
	typedef typename Graph::vertex_property_type type;
};

template <class G>
typename G::vertex_descriptor
add_vertex(const typename G::vertex_property_type& vp, G& g)
{
	return g.add_vertex(vp);
}

// EdgeMutablePropertyGraph

template <class Graph>
struct edge_property {
	typedef typename Graph::edge_property_type type;
};

template <class G>
std::pair<typename G::edge_descriptor, bool>
add_edge(
	typename G::vertex_descriptor u,
	typename G::vertex_descriptor v,
	const typename G::edge_property_type& ep,
	G& g)
{
	return g.add_edge(u, v, ep);
}

// PropertyGraph

template <class G>
typename vertex_property<G>::type
get(vertex_bundle_t, const G& g, typename G::vertex_descriptor u)
{
	return g[u];
}

template <class G>
typename edge_property<G>::type
get(edge_bundle_t, const G& g, typename G::edge_descriptor e)
{
	return g[e];
}

#include <istream>
#include <ostream>

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

#endif
