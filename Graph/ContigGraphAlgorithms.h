#ifndef CONTIGGRAPHALGORITHMS_H
#define CONTIGGRAPHALGORITHMS_H 1

#include "Graph.h"
#include <cassert>
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

/** Assemble an unambiguous path starting at vertex v. */
template<typename Graph, typename OutIt>
OutIt assemble(const Graph& g,
		typename Graph::vertex_descriptor v, OutIt out)
{
	for (; contiguous_out(g, v); v = *g.adjacent_vertices(v).first)
		*out++ = v;
	*out++ = v;
	return out;
}

#include "ContigGraph.h"
#include "ContigProperties.h"
#include <ostream>

void assemble(ContigGraph<ContigProperties>& g, std::ostream& out);

#endif
