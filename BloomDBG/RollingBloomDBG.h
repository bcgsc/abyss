/**
 * de Bruijn Graph data structure using a Bloom filter
 * Copyright 2015 Shaun Jackman, Ben Vandervalk
 */

#ifndef ROLLING_BLOOM_DBG_H
#define ROLLING_BLOOM_DBG_H 1

#include "Assembly/SeqExt.h" // for NUM_BASES
#include "Common/Hash.h"
#include "BloomDBG/MaskedKmer.h"
#include "Graph/Properties.h"
#include "BloomDBG/RollingHash.h"
#include "lib/bloomfilter-2dfba08d120d7659e8c75cf5c501b3b9040e98cb/BloomFilter.hpp"

#include <algorithm>
#include <cassert>
#include <cstdlib> // for abort
#include <fstream>
#include <string>
#include <utility> // for std::pair
#include <vector>
#include <iostream>

using boost::graph_traits;

NAMESPACE_STD_HASH_BEGIN
template <> struct hash< std::pair<MaskedKmer, RollingHash> > {
	/**
	 * Hash function for graph vertex type (vertex_descriptor),
	 * which is std::pair<MaskedKmer, RollingHash>.
	 */
	size_t operator()(const std::pair<MaskedKmer, RollingHash>& vertex) const
	{
		return vertex.second.getHash().at(0);
	}
};
NAMESPACE_STD_HASH_END

template <typename BF>
class RollingBloomDBG: public BF {
  public:
	/** The bundled vertex properties. */
	typedef no_property vertex_bundled;
	typedef no_property vertex_property_type;

	/** The bundled edge properties. */
	typedef no_property edge_bundled;
	typedef no_property edge_property_type;

	/** The bloom filter */
	const BF& m_bloom;

	RollingBloomDBG(const BF& bloom) : m_bloom(bloom) { }

  private:
	/** Copy constructor. */
	RollingBloomDBG(const RollingBloomDBG<BF>&);

}; // class RollingBloomDBG

// Graph

namespace boost {

/** Graph traits */
template <typename BF>
struct graph_traits< RollingBloomDBG<BF> > {
	// Graph

	/**
	 * Identifier for accessing a vertex in the graph.
	 * The second member of the pair (std::vector<size_t>) is
	 * a set of hash values associated with the k-mer.
	 */
	typedef std::pair<MaskedKmer, RollingHash> vertex_descriptor;
	typedef boost::directed_tag directed_category;
	struct traversal_category
		: boost::adjacency_graph_tag,
		boost::bidirectional_graph_tag,
		boost::vertex_list_graph_tag
		{ };
	typedef boost::disallow_parallel_edge_tag edge_parallel_category;

	static vertex_descriptor null_vertex()
	{
		return vertex_descriptor(MaskedKmer(), std::vector<size_t>());
	}

	// IncidenceGraph
	typedef std::pair<vertex_descriptor, vertex_descriptor>
		edge_descriptor;
	typedef unsigned degree_size_type;

	// VertexListGraph
	typedef size_t vertices_size_type;
	typedef void vertex_iterator;

	// EdgeListGraph
	typedef size_t edges_size_type;
	typedef void edge_iterator;

// AdjacencyGraph

/** Iterate through the adjacent vertices of a vertex. */
struct adjacency_iterator
	: public std::iterator<std::input_iterator_tag, vertex_descriptor>
{
	/** Skip to the next edge that is present. */
	void next()
	{
		for (; m_i < NUM_BASES; ++m_i) {
			m_v.first.setLastBase(SENSE, m_i);
			m_v.second = m_u.second;
			m_v.second.rollRight(
				m_u.first.getFirstBaseChar(),
				m_v.first.getLastBaseChar());
			if (vertex_exists(m_v, *m_g))
				break;
		}
	}

  public:

	adjacency_iterator() { }

	adjacency_iterator(const RollingBloomDBG<BF>& g) : m_g(&g), m_i(NUM_BASES) { }

	adjacency_iterator(const RollingBloomDBG<BF>& g, vertex_descriptor u)
		: m_g(&g), m_u(u), m_v(u), m_i(0)
	{
		m_v.first.shift(SENSE);
		next();
	}

	const vertex_descriptor& operator*() const
	{
		assert(m_i < NUM_BASES);
		return m_v;
	}

	bool operator==(const adjacency_iterator& it) const
	{
		return m_i == it.m_i;
	}

	bool operator!=(const adjacency_iterator& it) const
	{
		return !(*this == it);
	}

	adjacency_iterator& operator++()
	{
		assert(m_i < NUM_BASES);
		++m_i;
		next();
		return *this;
	}

	adjacency_iterator operator++(int)
	{
		adjacency_iterator it = *this;
		++*this;
		return it;
	}

  private:
	const RollingBloomDBG<BF>* m_g;
	vertex_descriptor m_u;
	vertex_descriptor m_v;
	short unsigned m_i;
}; // adjacency_iterator

/** IncidenceGraph */
struct out_edge_iterator
	: public std::iterator<std::input_iterator_tag, edge_descriptor>
{
	/** Skip to the next edge that is present. */
	void next()
	{
		for (; m_i < NUM_BASES; ++m_i) {
			m_v.first.setLastBase(SENSE, m_i);
			m_v.second = m_u.second;
			m_v.second.rollRight(
				m_u.first.getFirstBaseChar(),
				m_v.first.getLastBaseChar());
			if (vertex_exists(m_v, *m_g))
				break;
		}
	}

  public:
	out_edge_iterator() { }

	out_edge_iterator(const RollingBloomDBG<BF>& g) : m_g(&g), m_i(NUM_BASES) { }

	out_edge_iterator(const RollingBloomDBG<BF>& g, vertex_descriptor u)
		: m_g(&g), m_u(u), m_v(u), m_i(0)
	{
		m_v.first.shift(SENSE);
		next();
	}

	edge_descriptor operator*() const
	{
		assert(m_i < NUM_BASES);
		return edge_descriptor(m_u, m_v);
	}

	bool operator==(const out_edge_iterator& it) const
	{
		return m_i == it.m_i;
	}

	bool operator!=(const out_edge_iterator& it) const
	{
		return !(*this == it);
	}

	out_edge_iterator& operator++()
	{
		assert(m_i < NUM_BASES);
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
	const RollingBloomDBG<BF>* m_g;
	vertex_descriptor m_u;
	vertex_descriptor m_v;
	unsigned m_i;
}; // out_edge_iterator

/** BidirectionalGraph */
struct in_edge_iterator
	: public std::iterator<std::input_iterator_tag, edge_descriptor>
{
	/** Skip to the next edge that is present. */
	void next()
	{
		for (; m_i < NUM_BASES; ++m_i) {
			m_v.first.setLastBase(ANTISENSE, m_i);
			m_v.second = m_u.second;
			m_v.second.rollLeft(
				m_v.first.getFirstBaseChar(),
				m_u.first.getLastBaseChar());
			if (vertex_exists(m_v, *m_g))
				break;
		}
	}

  public:
	in_edge_iterator() { }

	in_edge_iterator(const RollingBloomDBG<BF>& g) : m_g(&g), m_i(NUM_BASES) { }

	in_edge_iterator(const RollingBloomDBG<BF>& g, vertex_descriptor u)
		: m_g(&g), m_u(u), m_v(u), m_i(0)
	{
		m_v.first.shift(ANTISENSE);
		next();
	}

	edge_descriptor operator*() const
	{
		assert(m_i < NUM_BASES);
		return edge_descriptor(m_v, m_u);
	}

	bool operator==(const in_edge_iterator& it) const
	{
		return m_i == it.m_i;
	}

	bool operator!=(const in_edge_iterator& it) const
	{
		return !(*this == it);
	}

	in_edge_iterator& operator++()
	{
		assert(m_i < NUM_BASES);
		++m_i;
		next();
		return *this;
	}

	in_edge_iterator operator++(int)
	{
		in_edge_iterator it = *this;
		++*this;
		return it;
	}

  private:
	const RollingBloomDBG<BF>* m_g;
	vertex_descriptor m_u;
	vertex_descriptor m_v;
	unsigned m_i;
}; // in_edge_iterator

}; // graph_traits<RollingBloomDBG>

} // namespace boost

// Subgraph

/** Return whether this vertex exists in the subgraph. */
template <typename Graph>
static inline bool
vertex_exists(typename graph_traits<Graph>::vertex_descriptor u, const Graph& g)
{
	return g.m_bloom.contains(u.second.getHash());
}

template <typename Graph>
static inline
std::pair<typename graph_traits<Graph>::adjacency_iterator,
		typename graph_traits<Graph>::adjacency_iterator>
adjacent_vertices(
		typename graph_traits<Graph>::vertex_descriptor u, const Graph& g)
{
	typedef typename graph_traits<Graph>::adjacency_iterator adjacency_iterator;
	return std::make_pair(adjacency_iterator(g, u), adjacency_iterator(g));
}

// IncidenceGraph
template <typename Graph>
static inline
typename graph_traits<Graph>::degree_size_type
out_degree(
		typename graph_traits<Graph>::vertex_descriptor u,
		const Graph& g)
{
	typedef typename graph_traits<Graph>::adjacency_iterator Ait;
	std::pair<Ait, Ait> adj = adjacent_vertices(u, g);
	return std::distance(adj.first, adj.second);
}

template <typename Graph>
static inline typename
std::pair<typename graph_traits<Graph>::out_edge_iterator,
	typename graph_traits<Graph>::out_edge_iterator>
out_edges(
		typename graph_traits<Graph>::vertex_descriptor u,
		const Graph& g)
{
	typedef typename graph_traits<Graph>::out_edge_iterator Oit;
	return std::make_pair(Oit(g, u), Oit(g));
}

// BidirectionalGraph
template <typename Graph>
static inline
std::pair<typename graph_traits<Graph>::in_edge_iterator,
	typename graph_traits<Graph>::in_edge_iterator>
in_edges(
		typename graph_traits<Graph>::vertex_descriptor u,
		const Graph& g)
{
	typedef typename graph_traits<Graph>::in_edge_iterator Iit;
	return std::make_pair(Iit(g, u), Iit(g));
}

template <typename Graph>
static inline
typename graph_traits<Graph>::degree_size_type
in_degree(typename graph_traits<Graph>::vertex_descriptor u,
		  const Graph& g)
{
	//return out_degree(reverseComplement(u), g);
	typedef typename graph_traits<Graph>::in_edge_iterator Iit;
	std::pair<Iit, Iit> it = in_edges(u, g);
	return std::distance(it.first, it.second);
}

// PropertyGraph

/** Return the reverse complement of the specified k-mer. */
template <typename Graph>
static inline
typename graph_traits<Graph>::vertex_descriptor
get(vertex_complement_t, const Graph&,
		typename graph_traits<Graph>::vertex_descriptor u)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	return V(reverseComplement(u.first), u.second);
}

/** Return the name of the specified vertex. */
template <typename Graph>
static inline
MaskedKmer get(vertex_name_t, const Graph&,
		typename graph_traits<Graph>::vertex_descriptor u)
{
	return u.first;
}

template <typename Graph>
static inline
bool
get(vertex_removed_t, const Graph&,
		typename graph_traits<Graph>::vertex_descriptor)
{
	return false;
}

template <typename Graph>
static inline
no_property
get(vertex_bundle_t, const Graph&,
		typename graph_traits<Graph>::edge_descriptor)
{
	return no_property();
}

template <typename Graph>
static inline
no_property
get(edge_bundle_t, const Graph&,
		typename graph_traits<Graph>::edge_descriptor)
{
	return no_property();
}

#endif
