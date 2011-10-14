/* Directed acyclic word graph. */
#ifndef DAWG_H
#define DAWG_H 1

#include "FMIndex.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <algorithm> // for distance
#include <cassert>
#include <utility>

using boost::graph_traits;

/** An FM-index is a representation of a DAWG. */
typedef FMIndex DAWG;

// Graph

namespace boost {

template <>
struct graph_traits<DAWG> {
	typedef DAWG::size_type size_type;

	// Graph
	typedef std::pair<size_type, size_type> vertex_descriptor;
	typedef boost::directed_tag directed_category;
	struct traversal_category : boost::incidence_graph_tag { };
	typedef boost::disallow_parallel_edge_tag edge_parallel_category;

	// IncidenceGraph
	typedef std::pair<vertex_descriptor, vertex_descriptor>
		edge_descriptor;
	typedef unsigned degree_size_type;

// VertexListGraph

/** Vertex iterator. */
struct vertex_iterator : public vertex_descriptor {
	vertex_iterator() : vertex_descriptor(0, 0) { }
	vertex_iterator(size_type l, size_type u)
		: vertex_descriptor(l, u) { }
	vertex_descriptor operator*() const { return *this; };
};

// IncidenceGraph

/** Out edge iterator. */
class out_edge_iterator
	: public std::iterator<std::input_iterator_tag, edge_descriptor>
{
  public:
	out_edge_iterator() : m_g(NULL), m_u(), m_v(), m_i(0) { }

	out_edge_iterator(const DAWG& g, vertex_descriptor u,
			degree_size_type i)
		: m_g(&g), m_u(u), m_v(), m_i(i)
	{
		next();
	}

	edge_descriptor operator*() const
	{
		return edge_descriptor(m_u, m_v);
	}

	bool operator==(const out_edge_iterator& it) const
	{
		return m_g == it.m_g && m_u == it.m_u && m_i == it.m_i;
	}

	bool operator!=(const out_edge_iterator& it) const
	{
		return !(*this == it);
	}

	out_edge_iterator& operator++()
	{
		assert(m_i < m_g->alphabetSize());
		++m_i;
		next();
		return *this;
	}

	out_edge_iterator operator++(int)
	{
		out_edge_iterator it = *this;
		++*this;
		return it;
	}

  private:
	/** Skip to the next edge that is present. */
	void next()
	{
		for (; m_i < m_g->alphabetSize(); m_i++) {
			m_v = vertex_descriptor(
					m_g->update(m_u.first, m_i),
					m_g->update(m_u.second, m_i));
			if (m_v.first < m_v.second)
				break;
		}
	}

	const DAWG* m_g;
	vertex_descriptor m_u;
	vertex_descriptor m_v;
	degree_size_type m_i;
}; // out_edge_iterator

}; // graph_traits<DAWG>

} // namespace boost

// IncidenceGraph

static inline
std::pair<
	graph_traits<DAWG>::out_edge_iterator,
	graph_traits<DAWG>::out_edge_iterator>
out_edges(
		graph_traits<DAWG>::vertex_descriptor u,
		const DAWG& g)
{
	typedef graph_traits<DAWG>::out_edge_iterator Eit;
	return std::pair<Eit, Eit>(
			Eit(g, u, 0),
			Eit(g, u, g.alphabetSize()));
}

static inline
graph_traits<DAWG>::degree_size_type
out_degree(
		graph_traits<DAWG>::vertex_descriptor u,
		const DAWG& g)
{
	typedef graph_traits<DAWG>::out_edge_iterator Eit;
	std::pair<Eit, Eit> it = out_edges(u, g);
	return std::distance(it.first, it.second);
}

// VertexListGraph

static inline
std::pair<graph_traits<DAWG>::vertex_iterator,
	graph_traits<DAWG>::vertex_iterator>
vertices(const DAWG& g)
{
	typedef graph_traits<DAWG>::vertex_iterator Vit;
	return std::pair<Vit, Vit>(Vit(0, g.size() + 1), Vit());
}

// PropertyGraph

static inline
DAWG::value_type
get(boost::edge_name_t,
		const DAWG& g,
		graph_traits<DAWG>::edge_descriptor e)
{
	graph_traits<DAWG>::vertex_descriptor v = target(e, g);
	assert(v.first < v.second);
	return g.symbolAt(v.first);
}

#endif
