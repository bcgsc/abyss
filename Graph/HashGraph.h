#ifndef HASH_GRAPH_H
#define HASH_GRAPH_H

#include "Common/UnorderedMap.h"
#include "Common/UnorderedSet.h"
#include "Graph/Properties.h"
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <cassert>
#include <sstream>

#define	DEFAULT_VERTEX_COLOR "gray"

template <class VertexType>
class HashGraph
{
public:

	typedef HashGraph<VertexType> Graph;
	typedef boost::graph_traits<Graph> GraphTraits;
	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::vertex_property_type vertex_property_type;
	typedef typename GraphTraits::edge_descriptor edge_descriptor;
	typedef typename GraphTraits::edge_property_type edge_property_type;
	typedef typename GraphTraits::out_edge_iterator out_edge_iterator;
	typedef typename GraphTraits::adjacency_iterator adjacency_iterator;
	typedef typename GraphTraits::degree_size_type degree_size_type;
	typedef typename GraphTraits::vertex_iterator vertex_iterator;
	typedef typename GraphTraits::vertices_size_type vertices_size_type;

	typedef typename GraphTraits::VertexList VertexList;
	typedef typename GraphTraits::VertexListIterator VertexListIterator;
	typedef typename GraphTraits::EdgeList EdgeList;

protected:

	typedef unordered_map<vertex_descriptor, VertexList,
		hash<vertex_descriptor> > VertexMap;
	typedef std::pair<vertex_descriptor, VertexList> VertexMapEntry;
	typedef unordered_map<vertex_descriptor, std::string> VertexColorMap;

	VertexMap m_vertices;
	VertexColorMap m_vertexColor;
	size_t m_numEdges;

public:

	typedef typename VertexMap::const_iterator VertexMapIterator;

	HashGraph() : m_numEdges(0) {}

	size_t approxMemSize()
	{
		size_t pointer_size = sizeof(void *);
		size_t entry_bytes = m_numEdges * sizeof(vertex_descriptor);
		size_t filled_bucket_bytes = m_vertices.size() *
			(sizeof(typename VertexMap::value_type) +
			 3 * pointer_size);
		size_t empty_bucket_bytes = size_t((1.0 - m_vertices.load_factor()) *
			m_vertices.bucket_count() * pointer_size);
		return entry_bytes + filled_bucket_bytes + empty_bucket_bytes;
	}

	void set_vertex_color(const vertex_descriptor& v, const std::string& color)
	{
		m_vertexColor[v] = color;
	}

	std::string get_vertex_color(const vertex_descriptor& v) const
	{
		typename VertexColorMap::const_iterator i = m_vertexColor.find(v);
		if (i == m_vertexColor.end())
			return DEFAULT_VERTEX_COLOR;
		return i->second;
	}

	std::pair<VertexMapIterator, VertexMapIterator>
	get_vertex_map_entries() const
	{
		return std::pair<VertexMapIterator, VertexMapIterator>
			(m_vertices.begin(), m_vertices.end());
	}

	vertices_size_type num_vertices() const
	{
		return m_vertices.size();
	}

	std::pair<VertexListIterator, VertexListIterator>
	get_successors(const vertex_descriptor& v) const
	{
		typename VertexMap::const_iterator i = m_vertices.find(v);
		assert(i != m_vertices.end());
		VertexListIterator begin = i->second.begin();
		VertexListIterator end = i->second.end();
		return std::pair<VertexListIterator, VertexListIterator>
			(begin, end);
	}

	degree_size_type
	out_degree(const vertex_descriptor& v) const
	{
		typename VertexMap::const_iterator i = m_vertices.find(v);
		if (i == m_vertices.end())
			return 0;
		return i->second.size();
	}

	typename VertexMap::iterator
	init_vertex(const vertex_descriptor &v)
	{
		bool inserted;
		VertexList successors;
		successors.reserve(1);
		VertexMapEntry entry(v, successors);
		typename VertexMap::iterator i;
		boost::tie(i, inserted) = m_vertices.insert(entry);
		assert(inserted);
		return i;
	}

	std::pair<edge_descriptor, bool>
	add_edge(const vertex_descriptor& u, const vertex_descriptor& v)
	{
		bool edgeInserted = false;
		typename VertexMap::iterator i = m_vertices.find(u);
		if (i == m_vertices.end())
			i = init_vertex(u);
		VertexList& successors = i->second;
		if (std::find(successors.begin(), successors.end(), v) == successors.end()) {
			successors.push_back(v);
			m_numEdges++;
			edgeInserted = true;
		}
		i = m_vertices.find(v);
		if (i == m_vertices.end())
			init_vertex(v);
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

	// BidirectionalGraph
	typedef void in_edge_iterator;

	// VertexListGraph
	typedef std::vector<vertex_descriptor> VertexList;
	typedef typename
		std::vector<vertex_descriptor>::const_iterator
		VertexListIterator;
	typedef unsigned vertices_size_type;

	// EdgeListGraph
	typedef void edges_size_type;

	// PropertyGraph
	typedef std::string vertex_bundled;
	typedef std::string vertex_property_type;
	typedef no_property edge_bundled;
	typedef no_property edge_property_type;

	// IncidenceGraph
	typedef std::vector<edge_descriptor> EdgeList;
	typedef unsigned degree_size_type;

	class vertex_iterator
		: public std::iterator<std::input_iterator_tag,
			const vertex_descriptor>
	{
		protected:

			void init(bool end)
			{
				boost::tie(m_it, m_end) = m_g->get_vertex_map_entries();
				if (end)
					m_it = m_end;
			}

		public:

			vertex_iterator() {}

			vertex_iterator(const HashGraph<VertexType>& g, bool end) : m_g(&g)
			{
				init(end);
			}

			vertex_descriptor operator *() const
			{
				return m_it->first;
			}

			bool operator ==(const vertex_iterator& it) const
			{
				return m_it == it.m_it;
			}

			bool operator !=(const vertex_iterator& it) const
			{
				return m_it != it.m_it;
			}

			vertex_iterator& operator ++()
			{
				++m_it;
				return *this;
			}

			vertex_iterator operator++(int)
			{
				vertex_iterator it = *this;
				++*this;
				return it;
			}

		private:

			const HashGraph<VertexType>* m_g;
			typename HashGraph<VertexType>::VertexMapIterator m_it, m_end;

	};

	struct adjacency_iterator
		: public std::iterator<std::input_iterator_tag, vertex_descriptor>
	{
		protected:

			void init(bool end)
			{
				boost::tie(m_successor, m_successor_end) = m_g->get_successors(m_v);
				if (end)
					m_successor = m_successor_end;
			}

		public:

			adjacency_iterator() { }

			adjacency_iterator(const HashGraph<VertexType>& g, vertex_descriptor v, bool end)
				: m_g(&g), m_v(v)
			{
				init(end);
			}

			vertex_descriptor operator*() const
			{
				return *m_successor;
			}

			bool operator==(const adjacency_iterator& it) const
			{
				return m_successor == it.m_successor;
			}

			bool operator!=(const adjacency_iterator& it) const
			{
				return !(*this == it);
			}

			adjacency_iterator& operator++()
			{
				m_successor++;
				return *this;
			}

			adjacency_iterator operator++(int)
			{
				adjacency_iterator it = *this;
				++*this;
				return it;
			}

		private:

			const HashGraph<VertexType>* m_g;
			vertex_descriptor m_v;
			VertexListIterator m_successor;
			VertexListIterator m_successor_end;
	};

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
		VertexListIterator m_successor;
		VertexListIterator m_successor_end;

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

// AdjacencyGraph

template <class VertexType>
std::pair<
	typename HashGraph<VertexType>::adjacency_iterator,
	typename HashGraph<VertexType>::adjacency_iterator>
adjacent_vertices(
	typename HashGraph<VertexType>::vertex_descriptor v,
	const HashGraph<VertexType>& g)
{
	typedef typename HashGraph<VertexType>::adjacency_iterator adjacency_iterator;
	return std::pair<adjacency_iterator, adjacency_iterator>
		(adjacency_iterator(g, v, false), adjacency_iterator(g, v, true));
}

// VertexListGraph

template <class VertexType>
std::pair<
	typename HashGraph<VertexType>::vertex_iterator,
	typename HashGraph<VertexType>::vertex_iterator>
vertices(const HashGraph<VertexType>& g)
{
	typedef typename HashGraph<VertexType>::vertex_iterator vertex_iterator;
	return std::pair<vertex_iterator, vertex_iterator>
		(vertex_iterator(g, false), vertex_iterator(g, true));
}

template <class VertexType>
typename HashGraph<VertexType>::vertices_size_type
num_vertices(const HashGraph<VertexType>& g)
{
	return g.num_vertices();
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

// PropertyGraph
// (dummy functions so that HashGraph
// works with DotIO.h routines)

template <class VertexType>
bool get(vertex_removed_t, const HashGraph<VertexType>&,
		typename HashGraph<VertexType>::vertex_descriptor)
{
	return false;
}

template <class VertexType>
void put(vertex_removed_t, HashGraph<VertexType>&,
		typename HashGraph<VertexType>::vertex_descriptor,
		bool)
{
	return;
}

template <class VertexType>
typename HashGraph<VertexType>::vertex_descriptor
get(vertex_name_t, const HashGraph<VertexType>&,
		typename HashGraph<VertexType>::vertex_descriptor v)
{
	return v;
}

template <class VertexType>
std::string
get(vertex_bundle_t, const HashGraph<VertexType>& g,
		typename HashGraph<VertexType>::vertex_descriptor v)
{
	const std::string& color = g.get_vertex_color(v);
	std::ostringstream ss;
	ss << "color=" << color;
	assert(ss);
	return ss.str();
}

#endif
