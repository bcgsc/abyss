#ifndef MAPGRAPH_H
#define MAPGRAPH_H 1

#include "ContigNode.h"
#include "Graph.h"
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

	/** Iterate through the vertices of this graph. */
	struct vertex_iterator : MapGraph::const_iterator
	{
		vertex_iterator(const MapGraph::const_iterator& it)
			: MapGraph::const_iterator(it) { }
		const vertex_descriptor& operator*() const
		{
			return MapGraph::const_iterator::operator*().first;
		}
	};

};

template <>
struct vertex_property<MapGraph> {
	typedef no_property type;
};

// IncidenceGraph

unsigned out_degree(graph_traits<MapGraph>::vertex_descriptor u,
		const MapGraph& g)
{
	MapGraph::const_iterator it = g.find(u);
	assert(it != g.end());
	return it->second.size();
}

// BidirectionalGraph

unsigned in_degree(graph_traits<MapGraph>::vertex_descriptor u,
		const MapGraph& g)
{
	return out_degree(~u, g);
}

// AdjacencyGraph

std::pair<graph_traits<MapGraph>::adjacency_iterator,
	graph_traits<MapGraph>::adjacency_iterator>
adjacent_vertices(graph_traits<MapGraph>::vertex_descriptor u,
		const MapGraph& g)
{
	MapGraph::const_iterator it = g.find(u);
	assert(it != g.end());
	return make_pair(it->second.begin(), it->second.end());
}

// VertexListGraph

std::pair<graph_traits<MapGraph>::vertex_iterator,
	graph_traits<MapGraph>::vertex_iterator>
vertices(const MapGraph& g)
{
	return make_pair(g.begin(), g.end());
}

// non-standard

static inline bool
is_removed(graph_traits<MapGraph>::vertex_descriptor, const MapGraph&)
{
	return false;
}

#endif
