#ifndef MAPGRAPH_H
#define MAPGRAPH_H 1

#include "ContigNode.h"
#include <map>
#include <set>
#include <utility> // for pair

typedef std::map<ContigNode, std::set<ContigNode> > MapGraph;

template <>
struct graph_traits<MapGraph> {
	typedef ContigNode vertex_descriptor;
	typedef MapGraph::mapped_type::const_iterator adjacency_iterator;
	typedef std::pair<vertex_descriptor, vertex_descriptor>
		edge_descriptor;
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

static inline graph_traits<MapGraph>::vertex_descriptor
source(graph_traits<MapGraph>::edge_descriptor u)
{
	return u.first;
}

static inline graph_traits<MapGraph>::vertex_descriptor
target(graph_traits<MapGraph>::edge_descriptor u)
{
	return u.second;
}

#endif
