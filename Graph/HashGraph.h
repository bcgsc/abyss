#ifndef HASH_GRAPH_H
#define HASH_GRAPH_H

#include "Common/UnorderedMap.h"
#include "Common/UnorderedSet.h"
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <cassert>
#include <iostream>

template <class VertexType>
class HashGraph
{
public:

	typedef HashGraph<VertexType> Graph;
	typedef boost::graph_traits<Graph> GraphTraits;
	typedef typename GraphTraits::EdgeList EdgeList;

	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::edge_descriptor edge_descriptor;
	typedef typename GraphTraits::out_edge_iterator out_edge_iterator;
	typedef typename GraphTraits::in_edge_iterator in_edge_iterator;
	typedef typename GraphTraits::degree_size_type degree_size_type;

	typedef unordered_map<vertex_descriptor, EdgeList> EdgeMap;
	typedef std::pair<vertex_descriptor, EdgeList> EdgeMapEntry;
	typedef typename EdgeMap::iterator EdgeMapIterator;

protected:

	unordered_set<vertex_descriptor> m_vertices;
	EdgeMap m_outEdges;
	EdgeMap m_inEdges;

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

	degree_size_type
	out_degree(vertex_descriptor v) const
	{
		typename EdgeMap::const_iterator mapi = m_outEdges.find(v);
		if (mapi == m_outEdges.end())
			return 0;
		else
			return mapi->second.size();
	}

	std::pair<in_edge_iterator, in_edge_iterator>
	in_edges(vertex_descriptor v) const
	{
		in_edge_iterator ei, ei_end;
		typename EdgeMap::const_iterator mapi = m_inEdges.find(v);

		if (mapi != m_inEdges.end()) {
			ei = mapi->second.begin();
			ei_end = mapi->second.end();
			std::cerr << "ei->second: " << ei->second << "\n";
		}

		return std::pair<in_edge_iterator, in_edge_iterator>(ei, ei_end);
	}

	degree_size_type
	in_degree(vertex_descriptor v) const
	{
		typename EdgeMap::const_iterator mapi = m_inEdges.find(v);
		if (mapi == m_inEdges.end())
			return 0;
		else
			return mapi->second.size();
	}

	std::pair<edge_descriptor, bool>
	add_edge(vertex_descriptor u, vertex_descriptor v)
	{
		edge_descriptor e(u, v);

		EdgeList& outEdges = m_outEdges[u];
		bool foundOutEdge = false;
		for (unsigned i = 0; i < outEdges.size(); i++) {
			if (outEdges[i] == e) {
				foundOutEdge = true;
				break;
			}
		}

		EdgeList& inEdges = m_inEdges[v];
		bool foundInEdge = false;
		for (unsigned i = 0; i < inEdges.size(); i++) {
			if (inEdges[i] == e) {
				foundInEdge = true;
				break;
			}
		}
		
		assert(foundInEdge == foundOutEdge);

		if (foundOutEdge)
			return std::pair<edge_descriptor, bool>(e, false);

		outEdges.push_back(e);
		inEdges.push_back(e);
		m_vertices.insert(u);
		m_vertices.insert(v);

		return std::pair<edge_descriptor, bool>(e, true);
	}

};

namespace boost {

template <class VertexType>
struct graph_traits<HashGraph<VertexType> > {

	// Graph
	typedef VertexType vertex_descriptor;

	// IncidenceGraph
	typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;
	typedef std::vector<edge_descriptor> EdgeList;
	typedef typename EdgeList::const_iterator edge_iterator;
	typedef edge_iterator out_edge_iterator;
	typedef unsigned degree_size_type;

	// BidirectionalGraph
	typedef edge_iterator in_edge_iterator;

	// AdjacencyGraph
	typedef void adjacency_iterator;

	// VertexListGraph
	typedef void vertices_size_type;

	// EdgeListGraph
	typedef void edges_size_type;
	typedef void vertex_iterator;

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

};

}

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

template <class VertexType>
typename HashGraph<VertexType>::degree_size_type
out_degree(
	typename HashGraph<VertexType>::vertex_descriptor v,
	const HashGraph<VertexType>& g)
{
	return g.out_degree(v);
}

// BidirectionalGraph

template <class VertexType>
std::pair<
	typename HashGraph<VertexType>::in_edge_iterator,
	typename HashGraph<VertexType>::in_edge_iterator>
in_edges(
	typename HashGraph<VertexType>::vertex_descriptor v,
	const HashGraph<VertexType>& g)
{
	return g.in_edges(v);
}

template <class VertexType>
typename HashGraph<VertexType>::degree_size_type
in_degree(
	typename HashGraph<VertexType>::vertex_descriptor v,
	const HashGraph<VertexType>& g)
{
	return g.in_degree(v);
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
