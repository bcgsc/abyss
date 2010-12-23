#ifndef CONTIGGRAPHALGORITHMS_H
#define CONTIGGRAPHALGORITHMS_H 1

#include "Algorithms.h"
#include "ContigNode.h"
#include "Functional.h"
#include "Graph.h"
#include "Iterator.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <set>
#include <utility>

/** Return true if the edge e is a palindrome. */
template<typename Graph>
struct IsPalindrome : std::unary_function<
		typename graph_traits<Graph>::edge_descriptor, bool>
{
	IsPalindrome(const Graph& g) : m_g(g) { }
	bool operator()(
			typename graph_traits<Graph>::edge_descriptor e) const
	{
		return source(e, m_g) == ~target(e, m_g);
	}
  private:
	const Graph& m_g;
};

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

/** Assemble a path of unambigous out edges starting at vertex u.
 * u itself is not copied to out.
 */
template<typename Graph, typename OutIt>
OutIt extend(const Graph& g,
		typename Graph::vertex_descriptor u, OutIt out)
{
	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	std::set<vertex_descriptor> seen;
	while (out_degree(u, g) == 1 && seen.insert(u).second) {
		u = *adjacent_vertices(u, g).first;
		*out++ = u;
	}
	return out;
}

/** Assemble an unambiguous path starting at vertex u.
 * Every edge must satisfy the predicate. */
template<typename Graph, typename OutIt, typename Predicate>
OutIt assemble_if(const Graph& g,
		typename Graph::vertex_descriptor u, OutIt out,
		Predicate pred)
{
	typedef typename graph_traits<Graph>::edge_descriptor
		edge_descriptor;
	while (contiguous_out(g, u)) {
		edge_descriptor e = *out_edges(u, g).first;
		if (!pred(e))
			break;
		*out++ = u;
		u = target(e, g);
	}
	*out++ = u;
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

/** Add the vertex and edge propeties of the path [first, last). */
template<typename Graph, typename It>
typename vertex_property<Graph>::type addProp(const Graph& g,
		It first, It last)
{
	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	typedef typename vertex_property<Graph>::type
		vertex_property_type;
	assert(first != last);
	vertex_property_type vp = get(vertex_bundle, g, *first);
	for (It it = first + 1; it != last; ++it) {
		vertex_descriptor u = *(it - 1);
		vertex_descriptor v = *it;
		vp += get(edge_bundle, g, u, v);
		vp += get(vertex_bundle, g, v);
	}
	return vp;
}

/** Merge the vertices in the sequence [first, last).
 * Create a new vertex whose property is the sum of [first, last).
 * Remove the vertices [first, last).
 */
template<typename Graph, typename It>
void merge(Graph& g, It first, It last)
{
	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	vertex_descriptor u = add_vertex(addProp(g, first, last), g);
	copy_in_edges(g, *first, u);
	copy_out_edges(g, *(last - 1), u);
}

/** Assemble unambiguous paths. Write the paths to out.
 * Every edge must satisfy the predicate. */
template<typename Graph, typename OutIt, typename Predicate>
OutIt assemble_if(Graph& g, OutIt out, Predicate pred0)
{
	typedef typename Graph::vertex_descriptor vertex_descriptor;
	typedef typename Graph::vertex_iterator vertex_iterator;
	// pred(e) = !isPalindrome(e) && pred0(e)
	binary_compose<std::logical_and<bool>,
		std::unary_negate<IsPalindrome<Graph> >, Predicate>
		pred(compose2(std::logical_and<bool>(),
				std::not1(IsPalindrome<Graph>(g)), pred0));
	std::pair<vertex_iterator, vertex_iterator> uit = g.vertices();
	for (vertex_iterator u = uit.first; u != uit.second; ++u) {
		if (!contiguous_out(g, *u) || contiguous_in(g, *u)
				|| !pred(*out_edges(*u, g).first))
			continue;
		typename output_iterator_traits<OutIt>::value_type path;
		assemble_if(g, *u, back_inserter(path), pred);
		assert(path.size() >= 2);
		assert(path.front() != path.back());
		merge(g, path.begin(), path.end());
		remove_vertex_if(g, path.begin(), path.end(),
				not1(std::mem_fun_ref(&ContigNode::ambiguous)));
		*out++ = path;
	}
	return out;
}

/** Assemble unambiguous paths. Write the paths to out. */
template<typename Graph, typename OutIt>
OutIt assemble(Graph& g, OutIt out)
{
	typedef typename graph_traits<Graph>::edge_descriptor
		edge_descriptor;
	return assemble_if(g, out, True<edge_descriptor>());
}

#endif
