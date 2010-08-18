#ifndef MAPGRAPH_H
#define MAPGRAPH_H 1

#include "ContigNode.h"
#include <map>
#include <set>

typedef std::map<ContigNode, std::set<ContigNode> > MapGraph;

template <>
struct graph_traits<MapGraph> {
	typedef ContigNode vertex_descriptor;
	typedef MapGraph::mapped_type::const_iterator adjacency_iterator;

	/** An edge (a pair of vertices). */
	struct edge_descriptor {
		edge_descriptor(const ContigNode& t, const ContigNode& h)
			: t(t), h(h) { }

		bool operator <(const edge_descriptor& o) const
		{
			return t != o.t ? t < o.t : h < o.h;
		}

		ContigNode t, h;
	};
};

unsigned out_degree(graph_traits<MapGraph>::vertex_descriptor u,
		const MapGraph& g)
{
	MapGraph::const_iterator it = g.find(u);
	assert(it != g.end());
	return it->second.size();
}

unsigned in_degree(graph_traits<MapGraph>::vertex_descriptor u,
		const MapGraph& g)
{
	return out_degree(~u, g);
}

std::pair<graph_traits<MapGraph>::adjacency_iterator,
	graph_traits<MapGraph>::adjacency_iterator>
adjacent_vertices(graph_traits<MapGraph>::vertex_descriptor u,
		const MapGraph& g)
{
	MapGraph::const_iterator it = g.find(u);
	assert(it != g.end());
	return make_pair(it->second.begin(), it->second.end());
}

#endif
