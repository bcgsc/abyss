#ifndef GRAPHALGORITHMS_H
#define GRAPHALGORITHMS_H 1

#include "Graph.h"
#include <cassert>
#include <set>
#include <vector>

/**
 * Remove transitive edges from the specified graph.
 * Find and remove the subset of edges (u,w) in E for which there
 * exists a vertex v such that the edges (u,v) and (v,w) exist in E.
 * This algorithm is not a general-purpose transitive reduction
 * algorithm. It is able to find transitive edges with exactly one
 * intermediate vertex.
 * @return the number of transitive edges removed from g
 */
template <typename Graph>
unsigned remove_transitive_edges(Graph& g)
{
	typedef typename graph_traits<Graph>::adjacency_iterator
		adjacency_iterator;
	typedef typename graph_traits<Graph>::edge_descriptor
		edge_descriptor;
	typedef typename graph_traits<Graph>::out_edge_iterator
		out_edge_iterator;
	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	typedef typename graph_traits<Graph>::vertex_iterator
		vertex_iterator;

	// The set of transitive edges.
	typedef std::vector<edge_descriptor> Edges;
	Edges transitive;

	std::pair<vertex_iterator, vertex_iterator> urange = vertices(g);
	for (vertex_iterator uit = urange.first;
			uit != urange.second; ++uit) {
		vertex_descriptor u = *uit;
		if (get(vertex_removed, g, u))
			continue;

		// Compute the set of vertices reachable in two hops.
		std::set<vertex_descriptor> seen;
		std::pair<adjacency_iterator, adjacency_iterator>
			vrange = adjacent_vertices(u, g);
		for (adjacency_iterator vit = vrange.first;
				vit != vrange.second; ++vit) {
			vertex_descriptor v = *vit;
			assert(u != v); // no self loops
			std::pair<adjacency_iterator, adjacency_iterator>
				wrange = adjacent_vertices(v, g);
			for (adjacency_iterator wit = wrange.first;
					wit != wrange.second; ++wit) {
				vertex_descriptor w = *wit;
				assert(v != w); // no self loops
				seen.insert(w);
			}
		}

		std::pair<out_edge_iterator, out_edge_iterator>
			uvrange = out_edges(u, g);
		for (out_edge_iterator uvit = uvrange.first;
				uvit != uvrange.second; ++uvit) {
			edge_descriptor uv = *uvit;
			vertex_descriptor v = target(uv, g);
			if (seen.count(v) > 0) {
				// The edge (u,v) is transitive. Mark it for removal.
				transitive.push_back(uv);
			}
		}
	}

	for (typename Edges::const_iterator it = transitive.begin();
			it != transitive.end(); ++it)
		remove_edge(*it, g);
	return transitive.size();
}

#endif
