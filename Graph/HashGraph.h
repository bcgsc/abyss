#ifndef HASH_GRAPH_H
#define HASH_GRAPH_H

#include "Common/UnorderedMap.h"
#include "Common/UnorderedSet.h"
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <cassert>

template <class VertexType>
class HashGraph
{
public:

	// Graph
	typedef VertexType vertex_descriptor;

	// IncidenceGraph
	typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;
	typedef unsigned degree_size_type;

	// BidirectionalGraph
	typedef void in_edge_iterator;

	// VertexListGraph
	typedef void vertices_size_type;

	// EdgeListGraph
	typedef void edges_size_type;

	// PropertyGraph
	typedef void vertex_bundled;
	typedef void vertex_property_type;
	typedef void edge_bundled;
	typedef void edge_property_type;

	typedef boost::directed_tag directed_category;
	typedef boost::allow_parallel_edge_tag edge_parallel_category;
	struct traversal_category
		: boost::incidence_graph_tag,
		boost::adjacency_graph_tag,
		boost::vertex_list_graph_tag,
		boost::edge_list_graph_tag { };

	typedef std::vector<edge_descriptor> EdgeList;

	typedef typename EdgeList::const_iterator out_edge_iterator;

protected:

	typedef unordered_map<vertex_descriptor, EdgeList> EdgeMap;
	typedef std::pair<vertex_descriptor, EdgeList> EdgeMapEntry;
	typedef typename EdgeMap::iterator EdgeMapIterator;

	unordered_set<vertex_descriptor> m_vertices;
	EdgeMap m_outEdges;

public:

	std::pair<out_edge_iterator, out_edge_iterator>
	out_edges(vertex_descriptor v) const
	{
		out_edge_iterator ei, ei_end;
		typename EdgeMap::const_iterator mapi = m_outEdges.find(v);

		if (mapi != m_outEdges.end()) {
			ei = mapi->second.begin();
			ei_end = mapi->second.end();
		}

		return std::pair<out_edge_iterator, out_edge_iterator>(ei, ei_end);
	}

	std::pair<edge_descriptor, bool>
	add_edge(vertex_descriptor u, vertex_descriptor v)
	{
		edge_descriptor e(u, v);

		EdgeList& outEdges = m_outEdges[u];
		for (unsigned i = 0; i < outEdges.size(); i++) {
			if (outEdges[i] == e) 
				return std::pair<edge_descriptor, bool>(e, false);
		}

		outEdges.push_back(e);
		return std::pair<edge_descriptor, bool>(e, true);
	}

};

// IncidenceGraph

template <class VertexType>
std::pair<
	typename HashGraph<VertexType>::out_edge_iterator,
	typename HashGraph<VertexType>::out_edge_iterator>
out_edges(
	typename HashGraph<VertexType>::vertex_descriptor v,
	const HashGraph<VertexType>& g)
{
	return g.out_edges(v);
}

// MutableGraph

template <class VertexType>
std::pair<typename HashGraph<VertexType>::edge_descriptor, bool>
add_edge(
	typename HashGraph<VertexType>::vertex_descriptor u,
	typename HashGraph<VertexType>::vertex_descriptor v,
	HashGraph<VertexType>& g)
{
	return g.add_edge(u, v);
}

#endif
