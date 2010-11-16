#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "Graph.h"
#include <utility>

/** A contig graph is a directed graph with the property that
 * the edge (u,v) implies the existence of the edge (~v,~u).
 */
template <typename G>
class ContigGraph : public G {
  public:
	// Graph
	typedef typename graph_traits<G>::vertex_descriptor
		vertex_descriptor;
	typedef typename graph_traits<G>::directed_category
		directed_category;
	typedef typename graph_traits<G>::traversal_category
		traversal_category;
	typedef typename graph_traits<G>::edge_parallel_category
		edge_parallel_category;

	// IncidenceGraph
	typedef typename graph_traits<G>::edge_descriptor
		edge_descriptor;
	typedef typename graph_traits<G>::out_edge_iterator
		out_edge_iterator;
	typedef typename graph_traits<G>::degree_size_type
		degree_size_type;

	// BidirectionalGraph
	typedef typename graph_traits<G>::in_edge_iterator
		in_edge_iterator;

	// AdjacencyGraph
	typedef typename graph_traits<G>::adjacency_iterator
		adjacency_iterator;

	// VertexListGraph
	typedef typename graph_traits<G>::vertex_iterator
		vertex_iterator;
	typedef typename graph_traits<G>::vertices_size_type
		vertices_size_type;

	// EdgeListGraph
	typedef typename graph_traits<G>::edge_iterator
		edge_iterator;
	typedef typename graph_traits<G>::edges_size_type
		edges_size_type;

	// VertexMutablePropertyGraph
	typedef typename vertex_property<G>::type vertex_property_type;

	// EdgeMutablePropertyGraph
	typedef typename edge_property<G>::type edge_property_type;

  public:
	/** Construct an empty contig graph. */
	ContigGraph() { }

	/** Construct a contig graph with n vertices. The underlying
	 * directed graph has two vertices for each contig. */
	ContigGraph(vertices_size_type n) : G(2 * n) { }

	/** Return the in degree of vertex v. */
	degree_size_type in_degree(vertex_descriptor v) const
	{
		return G::out_degree(~v);
	}

	/** Remove all out edges from vertex u. */
	void clear_out_edges(vertex_descriptor u)
	{
		std::pair<adjacency_iterator, adjacency_iterator>
			adj = G::adjacent_vertices(u);
		for (adjacency_iterator v = adj.first; v != adj.second; ++v) {
			if (~*v == u) {
				// When ~v == u, removing (~v,~u), which is (u,~u),
				// would invalidate our iterator. This edge will be
				// removed by clear_out_edges.
			} else
				G::remove_edge(~*v, ~u);
		}
		G::clear_out_edges(u);
	}

	/** Remove all in edges from vertex v. */
	void clear_in_edges(vertex_descriptor v)
	{
		clear_out_edges(~v);
	}

	/** Remove all edges to and from vertex v. */
	void clear_vertex(vertex_descriptor v)
	{
		clear_out_edges(v);
		clear_in_edges(v);
	}

	/** Add a vertex to this graph. */
	vertex_descriptor add_vertex(
			const vertex_property_type& data = vertex_property_type())
	{
		vertex_descriptor v = G::add_vertex(data);
		G::add_vertex(data);
		return v;
	}

	/** Remove vertex v from this graph. It is assumed that there
	 * are no edges to or from vertex v. It is best to call
	 * clear_vertex before remove_vertex.
	 */
	void remove_vertex(vertex_descriptor v)
	{
		G::remove_vertex(v);
		G::remove_vertex(~v);
	}

	/** Add edge (u,v) to this graph. */
	std::pair<edge_descriptor, bool>
	add_edge(vertex_descriptor u, vertex_descriptor v)
	{
		std::pair<edge_descriptor, bool> e = G::add_edge(u, v);
		if (u != ~v)
			G::add_edge(~v, ~u);
		return e;
	}

	/** Add edge (u,v) to this graph. */
	std::pair<edge_descriptor, bool>
	add_edge(vertex_descriptor u, vertex_descriptor v,
			const edge_property_type& ep)
	{
		std::pair<edge_descriptor, bool> e = G::add_edge(u, v, ep);
		if (u != ~v)
			G::add_edge(~v, ~u, ep);
		return e;
	}

	/** Remove the edge (u,v) from this graph. */
	void remove_edge(vertex_descriptor u, vertex_descriptor v)
	{
		G::remove_edge(u, v);
		if (u != ~v)
			G::remove_edge(~v, ~u);
	}

	/** Remove the edge e from this graph. */
	void remove_edge(edge_descriptor e)
	{
		remove_edge(source(e), target(e));
	}

  private:
	ContigGraph(const ContigGraph&);
};

#endif
