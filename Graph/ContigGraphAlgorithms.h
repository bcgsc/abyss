#ifndef CONTIGGRAPHALGORITHMS_H
#define CONTIGGRAPHALGORITHMS_H 1

#include "Algorithms.h"
#include "ContigNode.h"
#include "Graph.h"
#include "Iterator.h"
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

/** Add the outgoing edges of vertex u to vertex uout. */
template<typename Graph>
void copy_out_edges(Graph &g,
		typename Graph::vertex_descriptor u,
		typename Graph::vertex_descriptor uout)
{
	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	typedef typename graph_traits<Graph>::out_edge_iterator
		out_edge_iterator;
	typedef typename edge_property<Graph>::type edge_property_type;
	assert(u != uout);
	std::pair<out_edge_iterator, out_edge_iterator>
		edges = g.out_edges(u);
	bool palindrome = false;
	edge_property_type palindrome_ep;
	for (out_edge_iterator e = edges.first; e != edges.second; ++e) {
		vertex_descriptor v = target(*e, g);
		if (~v == u) {
			// When ~v == u, adding the edge (~v,~u), which is (u,~u),
			// would invalidate our iterator. Add the edge after this
			// loop completes.
			palindrome = true;
			palindrome_ep = g[*e];
		} else
			g.add_edge(uout, v, g[*e]);
	}
	if (palindrome) {
		g.add_edge(uout, ~u, palindrome_ep);
		g.add_edge(uout, ~uout, palindrome_ep);
	}
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

/** Remove vertices in the sequence [first, last) from the graph
 * for which the predicate p is true. Edges incident to those vertices
 * are removed as well.
 */
template<typename Graph, typename It, typename Predicate>
void remove_vertex_if(Graph& g, It first, It last, Predicate p)
{
	for_each_if(first, last,
			bind1st(std::mem_fun(&Graph::clear_vertex), &g), p);
	for_each_if(first, last,
			bind1st(std::mem_fun(&Graph::remove_vertex), &g), p);
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
}

/** Assemble unambiguous paths. Write the paths to out. */
template<typename Graph, typename OutIt>
OutIt assemble(Graph& g, OutIt out)
{
	typedef typename Graph::vertex_descriptor vertex_descriptor;
	typedef typename Graph::vertex_iterator vertex_iterator;
	std::pair<vertex_iterator, vertex_iterator> vit = g.vertices();
	for (vertex_iterator v = vit.first; v != vit.second; ++v) {
		if (!contiguous_out(g, *v) || contiguous_in(g, *v)
				|| is_palindrome(g, *out_edges(*v, g).first))
			continue;
		typename output_iterator_traits<OutIt>::value_type path;
		assemble(g, *v, back_inserter(path));
		assert(path.size() >= 2);
		assert(path.front() != path.back());
		merge(g, path.begin(), path.end());
		remove_vertex_if(g, path.begin(), path.end(),
				not1(std::mem_fun_ref(&ContigNode::ambiguous)));
		*out++ = path;
	}
	return out;
}

#endif
