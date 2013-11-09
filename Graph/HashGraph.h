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
	typedef typename GraphTraits::VertexList VertexList;
	typedef typename GraphTraits::EdgeList EdgeList;

	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::edge_descriptor edge_descriptor;
	typedef typename GraphTraits::out_edge_iterator out_edge_iterator;
	typedef typename GraphTraits::degree_size_type degree_size_type;
	typedef typename GraphTraits::vertex_iterator vertex_iterator;

protected:

	typedef unordered_map<vertex_descriptor, VertexList> VertexMap;
	VertexMap m_vertices;

public:

	std::pair<vertex_iterator, vertex_iterator>
	get_successors(const vertex_descriptor& v) const
	{
		typename VertexMap::const_iterator i = m_vertices.find(v);
		vertex_iterator begin, end;
		if (i != m_vertices.end()) {
			begin = i->second.begin();
			end = i->second.end();
		}
		return std::pair<vertex_iterator, vertex_iterator> (begin, end);
	}

	degree_size_type
	out_degree(const vertex_descriptor& v) const
	{
		typename VertexMap::const_iterator i = m_vertices.find(v);
		if (i == m_vertices.end())
			return 0;
		return i->second.size();
	}

	std::pair<edge_descriptor, bool>
	add_edge(const vertex_descriptor& u, const vertex_descriptor& v)
	{
		bool edgeInserted = false;
		VertexList& successors = m_vertices[u];
		if (std::find(successors.begin(), successors.end(), v) == successors.end()) {
			successors.push_back(v);
			edgeInserted = true;
		}
		return std::pair<edge_descriptor, bool>
			(edge_descriptor(u, v), edgeInserted);
	}

};

namespace boost {

template <class VertexType>
struct graph_traits< HashGraph<VertexType> > {

	// Graph
	typedef VertexType vertex_descriptor;
	typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;
	typedef boost::directed_tag directed_category;
	typedef boost::allow_parallel_edge_tag edge_parallel_category;
	struct traversal_category
		: boost::incidence_graph_tag,
		boost::adjacency_graph_tag,
		boost::vertex_list_graph_tag,
		boost::edge_list_graph_tag { };

	// AdjacencyGraph
	typedef void adjacency_iterator;

	// BidirectionalGraph
	typedef void in_edge_iterator;

	// VertexListGraph
	typedef std::vector<vertex_descriptor> VertexList;
	typedef typename VertexList::const_iterator vertex_iterator;
	typedef unsigned vertices_size_type;

	// EdgeListGraph
	typedef void edges_size_type;

	// PropertyGraph
	typedef void vertex_bundled;
	typedef void vertex_property_type;
	typedef void edge_bundled;
	typedef void edge_property_type;

	// IncidenceGraph
	typedef std::vector<edge_descriptor> EdgeList;
	typedef unsigned degree_size_type;

	struct out_edge_iterator
		: public std::iterator<std::input_iterator_tag, edge_descriptor>
	{

	protected:

		void init(bool end)
		{
			boost::tie(m_successor, m_successor_end) = m_g->get_successors(m_v);
			if (end)
				m_successor = m_successor_end;
		}

	public:

		out_edge_iterator() { }

		out_edge_iterator(const HashGraph<VertexType>& g, vertex_descriptor v)
			: m_g(&g), m_v(v)
		{
			init(false);
		}

		out_edge_iterator(const HashGraph<VertexType>& g, vertex_descriptor v, bool end)
			: m_g(&g), m_v(v)
		{
			init(end);
		}

		edge_descriptor operator*() const
		{
			return edge_descriptor(m_v, *m_successor);
		}

		bool operator==(const out_edge_iterator& it) const
		{
			return m_successor == it.m_successor;
		}

		bool operator!=(const out_edge_iterator& it) const
		{
			return !(*this == it);
		}

		out_edge_iterator& operator++()
		{
			m_successor++;
			return *this;
		}

		out_edge_iterator operator++(int)
		{
			out_edge_iterator it = *this;
			++*this;
			return it;
		}

	private:

		const HashGraph<VertexType>* m_g;
		vertex_descriptor m_v;
		vertex_iterator m_successor;
		vertex_iterator m_successor_end;

	}; // out_edge_iterator

}; // graph_traits

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
	typedef typename HashGraph<VertexType>::out_edge_iterator out_edge_iterator;
	return std::pair<out_edge_iterator, out_edge_iterator>
		(out_edge_iterator(g, v), out_edge_iterator(g, v, true));
}

template <class VertexType>
typename HashGraph<VertexType>::degree_size_type
out_degree(
	typename HashGraph<VertexType>::vertex_descriptor v,
	const HashGraph<VertexType>& g)
{
	return g.out_degree(v);
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
