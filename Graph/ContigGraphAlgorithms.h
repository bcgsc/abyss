#ifndef CONTIGGRAPHALGORITHMS_H
#define CONTIGGRAPHALGORITHMS_H 1

#include "Common/Algorithms.h"
#include "Common/ContigNode.h"
#include "Common/Estimate.h" // for BetterDistanceEst
#include "Common/Functional.h"
#include "Common/Iterator.h"
#include "Graph/ContigGraph.h"
#include <algorithm>
#include <boost/graph/graph_traits.hpp>
#include <cassert>
#include <functional>
#include <set>
#include <utility>

using boost::graph_traits;

/** Return true if the edge e is a palindrome. */
template<typename Graph>
struct IsPalindrome : std::unary_function<typename graph_traits<Graph>::edge_descriptor, bool>
{
	IsPalindrome(const Graph& g)
	  : m_g(g)
	{}
	bool operator()(typename graph_traits<Graph>::edge_descriptor e) const
	{
		return source(e, m_g) == get(vertex_complement, m_g, target(e, m_g));
	}

  private:
	const Graph& m_g;
};

/** Return whether the outgoing edge of vertex u is contiguous. */
template<typename Graph>
bool
contiguous_out(const Graph& g, typename graph_traits<Graph>::vertex_descriptor u)
{
	return out_degree(u, g) == 1 && in_degree(*adjacent_vertices(u, g).first, g) == 1;
}

/** Return whether the incoming edge of vertex u is contiguous. */
template<typename Graph>
bool
contiguous_in(const Graph& g, typename graph_traits<Graph>::vertex_descriptor u)
{
	return contiguous_out(g, get(vertex_complement, g, u));
}

/** Add the outgoing edges of vertex u to vertex uout. */
template<typename Graph>
void
copy_out_edges(
    Graph& g,
    typename Graph::vertex_descriptor u,
    typename Graph::vertex_descriptor uout)
{
	typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;
	typedef typename edge_property<Graph>::type edge_property_type;
	assert(u != uout);
	std::pair<out_edge_iterator, out_edge_iterator> edges = g.out_edges(u);
	bool palindrome = false;
	edge_property_type palindrome_ep;
	for (out_edge_iterator e = edges.first; e != edges.second; ++e) {
		vertex_descriptor v = target(*e, g);
		vertex_descriptor vc = get(vertex_complement, g, v);
		if (vc == u) {
			// When ~v == u, adding the edge (~v,~u), which is (u,~u),
			// would invalidate our iterator. Add the edge after this
			// loop completes.
			palindrome = true;
			palindrome_ep = g[*e];
		} else
			g.add_edge(uout, v, g[*e]);
	}
	if (palindrome) {
		vertex_descriptor uc = get(vertex_complement, g, u);
		vertex_descriptor uoutc = get(vertex_complement, g, uout);
		g.add_edge(uout, uc, palindrome_ep);
		g.add_edge(uout, uoutc, palindrome_ep);
	}
}

/** Add the incoming edges of vertex u to vertex v. */
template<typename Graph>
void
copy_in_edges(Graph& g, typename Graph::vertex_descriptor u, typename Graph::vertex_descriptor v)
{
	copy_out_edges(g, get(vertex_complement, g, u), get(vertex_complement, g, v));
}

/** Assemble a path of unambigous out edges starting at vertex u.
 * u itself is not copied to out.
 */
template<typename Graph, typename OutIt>
OutIt
extend(const Graph& g, typename Graph::vertex_descriptor u, OutIt out)
{
	typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
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
OutIt
assemble_if(const Graph& g, typename Graph::vertex_descriptor u, OutIt out, Predicate pred)
{
	typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
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
void
remove_vertex_if(Graph& g, It first, It last, Predicate p)
{
	for_each_if(
	    first, last, [&g](const ContigNode& c) { return g.clear_vertex(c); }, p);
	for_each_if(
	    first, last, [&g](const ContigNode& c) { return g.remove_vertex(c); }, p);
}

/** Add the vertex and edge propeties of the path [first, last). */
template<typename Graph, typename It, typename VP>
VP
addProp(const Graph& g, It first, It last, const VP*)
{
	typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	assert(first != last);
	VP vp = get(vertex_bundle, g, *first);
	for (It it = first + 1; it != last; ++it) {
		vertex_descriptor u = *(it - 1);
		vertex_descriptor v = *it;
		vp += get(edge_bundle, g, u, v);
		vp += get(vertex_bundle, g, v);
	}
	return vp;
}

template<typename Graph, typename It>
no_property
addProp(const Graph&, It, It, const no_property*)
{
	return no_property();
}

template<typename Graph, typename It>
typename vertex_property<Graph>::type
addProp(const Graph& g, It first, It last)
{
	return addProp(g, first, last, (typename vertex_property<Graph>::type*)NULL);
}

/** Merge the vertices in the sequence [first, last).
 * Create a new vertex whose property is the sum of [first, last).
 * Remove the vertices [first, last).
 */
template<typename Graph, typename It>
typename graph_traits<Graph>::vertex_descriptor
merge(Graph& g, It first, It last)
{
	typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	assert(first != last);
	vertex_descriptor u = add_vertex(addProp(g, first, last), g);
	copy_in_edges(g, *first, u);
	copy_out_edges(g, *(last - 1), u);
	return u;
}

/** Assemble unambiguous paths. Write the paths to out.
 * Every edge must satisfy the predicate. */
template<typename Graph, typename OutIt, typename Predicate>
OutIt
assemble_if(Graph& g, OutIt out, Predicate pred0)
{
	typedef typename Graph::vertex_iterator vertex_iterator;
	// pred(e) = !isPalindrome(e) && pred0(e)
	binary_compose<std::logical_and<bool>, std::unary_negate<IsPalindrome<Graph>>, Predicate> pred(
	    compose2(std::logical_and<bool>(), std::not1(IsPalindrome<Graph>(g)), pred0));
	std::pair<vertex_iterator, vertex_iterator> uit = g.vertices();
	for (vertex_iterator u = uit.first; u != uit.second; ++u) {
		if (!contiguous_out(g, *u) || contiguous_in(g, *u) || !pred(*out_edges(*u, g).first))
			continue;
		typename output_iterator_traits<OutIt>::value_type path;
		assemble_if(g, *u, back_inserter(path), pred);
		assert(path.size() >= 2);
		assert(path.front() != path.back());
		merge(g, path.begin(), path.end());
		remove_vertex_if(
		    g, path.begin(), path.end(), [](const ContigNode& c) { return !c.ambiguous(); });
		*out++ = path;
	}
	return out;
}

/** Assemble unambiguous paths. Write the paths to out. */
template<typename Graph, typename OutIt>
OutIt
assemble(Graph& g, OutIt out)
{
	typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
	return assemble_if(g, out, True<edge_descriptor>());
}

/** Return true if the edge e is +ve sense. */
template<typename Graph>
struct IsPositive : std::unary_function<typename graph_traits<Graph>::edge_descriptor, bool>
{
	IsPositive(const Graph& g)
	  : m_g(g)
	{}
	bool operator()(typename graph_traits<Graph>::edge_descriptor e) const
	{
		return !get(vertex_sense, m_g, source(e, m_g)) && !get(vertex_sense, m_g, target(e, m_g));
	}

  private:
	const Graph& m_g;
};

/** Assemble unambiguous paths in forward orientation only.
 * Write the paths to out. */
template<typename Graph, typename OutIt>
OutIt
assemble_stranded(Graph& g, OutIt out)
{
	return assemble_if(g, out, IsPositive<Graph>(g));
}

/** Remove tips.
 * For an edge (u,v), remove the vertex v if deg+(u) > 1,
 * deg+(v) = 0, and p(v) is true.
 * Stores all removed vertices in result.
 */
template<typename Graph, typename OutputIt, typename Pred>
OutputIt
pruneTips_if(Graph& g, OutputIt result, Pred p)
{
	typedef typename graph_traits<Graph>::adjacency_iterator Vit;
	typedef typename graph_traits<Graph>::vertex_iterator Uit;
	typedef typename graph_traits<Graph>::vertex_descriptor V;

	/** Identify the tips. */
	std::vector<V> tips;
	std::pair<Uit, Uit> urange = vertices(g);
	for (Uit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		if (out_degree(u, g) < 2)
			continue;
		std::pair<Vit, Vit> vrange = adjacent_vertices(u, g);
		for (Vit vit = vrange.first; vit != vrange.second; ++vit) {
			V v = *vit;
			// assert(v != u);
			if (out_degree(v, g) == 0 && p(v))
				tips.push_back(v);
		}
	}

	/** Remove the tips. */
	remove_vertex_if(g, tips.begin(), tips.end(), True<V>());
	std::transform(
	    tips.begin(), tips.end(), result, [](const ContigNode& c) { return c.contigIndex(); });
	return result;
}

/** Return true if the vertex is a normal 1-in 0-out tip. */
template<typename Graph>
struct IsTip : std::unary_function<typename graph_traits<Graph>::vertex_descriptor, bool>
{
	IsTip(const Graph& g)
	  : m_g(g)
	{}
	bool operator()(typename graph_traits<Graph>::vertex_descriptor v) const
	{
		return in_degree(v, m_g) == 1;
	}

  private:
	const Graph& m_g;
};

/** Remove tips.
 * For an edge (u,v), remove the vertex v if deg+(u) > 1
 * and deg-(v) = 1 and deg+(v) = 0.
 * Stores all removed vertices in result.
 */
template<typename Graph, typename OutputIt>
OutputIt
pruneTips(Graph& g, OutputIt result)
{
	return pruneTips_if(g, result, IsTip<Graph>(g));
}

/** Remove islands.
 * For a vertex v, remove v if deg+(v) = 0, deg-(v) = 0 and p(v) is
 * true.
 * Stores all removed vertices in result.
 */
template<typename Graph, typename OutputIt, typename Pred>
OutputIt
removeIslands_if(Graph& g, OutputIt result, Pred p)
{
	typedef typename graph_traits<Graph>::vertex_iterator Uit;
	typedef typename graph_traits<Graph>::vertex_descriptor V;

	/** Identify and remove Islands. */
	std::pair<Uit, Uit> urange = vertices(g);
	for (Uit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		if (get(vertex_removed, g, u))
			continue;
		if (p(u) && in_degree(u, g) == 0 && out_degree(u, g) == 0) {
			*result++ = get(vertex_contig_index, g, u);
			clear_vertex(u, g);
			remove_vertex(u, g);
		}
	}
	return result;
}

/** Add missing complementary edges. */
template<typename DG>
size_t
addComplementaryEdges(ContigGraph<DG>& g)
{
	typedef ContigGraph<DG> Graph;
	typedef graph_traits<Graph> GTraits;
	typedef typename GTraits::edge_descriptor E;
	typedef typename GTraits::edge_iterator Eit;
	typedef typename GTraits::vertex_descriptor V;

	std::pair<Eit, Eit> erange = edges(g);
	size_t numAdded = 0;
	for (Eit eit = erange.first; eit != erange.second; ++eit) {
		E e = *eit;
		V u = source(e, g), v = target(e, g);
		V uc = get(vertex_complement, g, u);
		V vc = get(vertex_complement, g, v);
		E f;
		bool found;
		boost::tie(f, found) = edge(vc, uc, g);
		if (!found) {
			add_edge(vc, uc, g[e], static_cast<DG&>(g));
			numAdded++;
		} else if (!(g[e] == g[f])) {
			// The edge properties do not agree. Select the better.
			g[e] = g[f] = BetterDistanceEst()(g[e], g[f]);
		}
	}
	return numAdded;
}

#endif
