#ifndef CONTIGGRAPHALGORITHMS_H
#define CONTIGGRAPHALGORITHMS_H 1

#include "Algorithms.h"
#include "ContigNode.h"
#include "Graph.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>
#include <utility>

/** Return whether the outgoing edge of vertex u is contiguous. */
template<typename Graph>
bool contiguous_out(const Graph& g,
		typename graph_traits<Graph>::vertex_descriptor u)
{
	return out_degree(u, g) == 1
		&& in_degree(*adjacent_vertices(u, g).first, g) == 1;
}

/** Return whether the incoming edge of vertex u is contiguous. */
template<typename Graph>
bool contiguous_in(const Graph& g,
		typename Graph::vertex_descriptor u)
{
	return contiguous_out(g, ~u);
}

/** Add the outgoing edges of vertex u to vertex v. */
template<typename Graph>
void copy_out_edges(Graph &g,
		typename Graph::vertex_descriptor u,
		typename Graph::vertex_descriptor v)
{
	typedef typename Graph::adjacency_iterator adjacency_iterator;
	assert(u != v);
	std::pair<adjacency_iterator, adjacency_iterator>
		adj = g.adjacent_vertices(u);
	for (adjacency_iterator it = adj.first; it != adj.second; ++it)
		g.add_edge(v, *it);
}

/** Add the incoming edges of vertex u to vertex v. */
template<typename Graph>
void copy_in_edges(Graph& g,
		typename Graph::vertex_descriptor u,
		typename Graph::vertex_descriptor v)
{
	copy_out_edges(g, ~u, ~v);
}

/** Return true if the edge e is palindromic. */
template<typename Graph>
bool is_palindrome(const Graph& g, typename Graph::edge_descriptor e)
{
	return source(e, g) == ~target(e, g);
}

/** Assemble an unambiguous path starting at vertex v. */
template<typename Graph, typename OutIt>
OutIt assemble(const Graph& g,
		typename Graph::vertex_descriptor v, OutIt out)
{
	for (; contiguous_out(g, v)
			&& !is_palindrome(g, *out_edges(v, g).first);
			v = *g.adjacent_vertices(v).first)
		*out++ = v;
	*out++ = v;
	return out;
}

template<typename Graph>
struct AddVertexProp {
	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	typedef typename vertex_property<Graph>::type
		vertex_property_type;
	const Graph& g;
	AddVertexProp(const Graph& g) : g(g) { }
	vertex_property_type operator()(
			const vertex_property_type& vp,
			const vertex_descriptor& u) const
	{
		return vp + get(vertex_bundle, g, u);
	}
};

/** Merge the vertices in the sequence [first, last).
 * Create a new vertex whose property is the sum of [first, last).
 * Remove the vertices [first, last).
 */
template<typename Graph, typename It>
void merge(Graph& g, It first, It last)
{
	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	typedef typename vertex_property<Graph>::type
		vertex_property_type;

	assert(first != last);
	vertex_property_type vp = std::accumulate(first + 1, last,
			get(vertex_bundle, g, *first), AddVertexProp<Graph>(g));
	vertex_descriptor u = add_vertex(vp, g);
	copy_in_edges(g, *first, u);
	copy_out_edges(g, *(last - 1), u);
	for_each_if(first, last,
			bind1st(std::mem_fun(&Graph::clear_vertex), &g),
			not1(std::mem_fun_ref(&ContigNode::ambiguous)));
	for_each_if(first, last,
			bind1st(std::mem_fun(&Graph::remove_vertex), &g),
			not1(std::mem_fun_ref(&ContigNode::ambiguous)));
}

#include "ContigGraph.h"
#include "ContigProperties.h"
#include "DirectedGraph.h"
#include <ostream>

void assemble(ContigGraph<DirectedGraph<ContigProperties> >& g,
		std::ostream& out);

#endif
