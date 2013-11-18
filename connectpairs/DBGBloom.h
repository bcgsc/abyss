/**
 * de Bruijn Graph data structure using a Bloom filter
 * Copyright 2013 Shaun Jackman
 */

#ifndef DBGBLOOM_H
#define DBGBLOOM_H 1

#include "CountingBloomFilter.h"

#include "Common/IOUtil.h"
#include "Common/Kmer.h"
#include "Common/SeqExt.h" // for NUM_BASES
#include "Graph/Properties.h"

#include <algorithm>
#include <cassert>
#include <cstdlib> // for abort
#include <fstream>
#include <string>

using boost::graph_traits;

/** de Bruijn Graph data structure using a Bloom filter. */
class DBGBloom {
  public:
	/** The bundled vertex properties. */
	typedef no_property vertex_bundled;
	typedef no_property vertex_property_type;

	/** The bundled edge properties. */
	typedef no_property edge_bundled;
	typedef no_property edge_property_type;

	/** The size of a k-mer */
	unsigned m_k;

	/** The contents of the FASTA file. */
	std::string m_fa;

	/** The bloom filter */
	CountingBloomFilter m_bloom;

	/** Construct a graph.
	 * @param k the size of a k-mer
	 * @param n the size of the bloom filter in bits
	 */
	DBGBloom(unsigned k, size_t n)
		: m_k(k), m_bloom(n) { }

	/** Load the Bloom filter from a string. */
	void assign(const std::string& s)
	{
		m_fa = s;
		for (size_t i = 0; i < s.size() - m_k + 1; ++i) {
			std::string kmer = s.substr(i, m_k);
			size_t pos = kmer.find_last_not_of("ACGTacgt");
			if (pos == std::string::npos) {
				Kmer u(kmer);
				m_bloom.insert(u);
				u.reverseComplement();
				m_bloom.insert(u);
			} else
				i += pos;
		}
	}

	/** Load the Bloom filter from a file. */
	void open(const std::string& path)
	{
		assert(!path.empty());
		readFile(path.c_str(), m_fa);
		assign(m_fa);
	}

  private:
	/** Copy constructor. */
	DBGBloom(const DBGBloom&);
}; // class DBGBloom

/** PropertyGraph vertex_index */
struct DBGBloomIndexMap
	: boost::put_get_helper<size_t, DBGBloomIndexMap>
{
	typedef Kmer key_type;
	typedef size_t value_type;
	typedef value_type reference;
	typedef boost::readable_property_map_tag category;

	const DBGBloom& m_g;

	DBGBloomIndexMap(const DBGBloom& g) : m_g(g) { }

	reference operator[](const key_type& u) const
	{
		return m_g.m_bloom.hash(u) % m_g.m_bloom.size();
	}
};

// Graph

namespace boost {

/** Graph traits */
template <>
struct graph_traits<DBGBloom> {
	// Graph
	typedef Kmer vertex_descriptor;
	typedef boost::directed_tag directed_category;
	struct traversal_category
		: boost::adjacency_graph_tag,
		boost::bidirectional_graph_tag,
		boost::vertex_list_graph_tag
		{ };
	typedef boost::disallow_parallel_edge_tag edge_parallel_category;

	static vertex_descriptor null_vertex() { return Kmer(); }

	// IncidenceGraph
	typedef std::pair<vertex_descriptor, vertex_descriptor>
		edge_descriptor;
	typedef unsigned degree_size_type;

	// VertexListGraph
	typedef size_t vertices_size_type;

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
			m_v.setLastBase(SENSE, m_i);
			if (vertex_exists(m_v, m_g))
				break;
		}
	}

  public:
	adjacency_iterator(const DBGBloom& g) : m_g(g), m_i(NUM_BASES) { }

	adjacency_iterator(const DBGBloom& g, vertex_descriptor u)
		: m_g(g), m_v(u), m_i(0)
	{
		m_v.shift(SENSE);
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

  private:
	const DBGBloom& m_g;
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
			m_v.setLastBase(SENSE, m_i);
			if (vertex_exists(m_v, *m_g))
				break;
		}
	}

  public:
	out_edge_iterator() { }

	out_edge_iterator(const DBGBloom& g) : m_g(&g), m_i(NUM_BASES) { }

	out_edge_iterator(const DBGBloom& g, vertex_descriptor u)
		: m_g(&g), m_u(u), m_v(u), m_i(0)
	{
		m_v.shift(SENSE);
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
	const DBGBloom* m_g;
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
			m_v.setLastBase(ANTISENSE, m_i);
			if (vertex_exists(m_v, *m_g))
				break;
		}
	}

  public:
	in_edge_iterator() { }

	in_edge_iterator(const DBGBloom& g) : m_g(&g), m_i(NUM_BASES) { }

	in_edge_iterator(const DBGBloom& g, vertex_descriptor u)
		: m_g(&g), m_u(u), m_v(u), m_i(0)
	{
		m_v.shift(ANTISENSE);
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
	const DBGBloom* m_g;
	vertex_descriptor m_u;
	vertex_descriptor m_v;
	unsigned m_i;
}; // in_edge_iterator

// Subgraph

/** Return whether this vertex exists in the subgraph. */
static inline
bool
vertex_exists(graph_traits<DBGBloom>::vertex_descriptor u,
		const DBGBloom& g)
{
	return g.m_bloom[u];
}

// VertexListGraph

/** Iterate through the vertices of this graph. */
struct vertex_iterator
	: public std::iterator<std::input_iterator_tag, vertex_descriptor>
{
	/** Skip to the next vertex that is present. */
	void next()
	{
		for (; m_i < m_g->m_fa.size() - m_g->m_k + 1; ++m_i) {
			std::string kmer = m_g->m_fa.substr(m_i, m_g->m_k);
			size_t pos = kmer.find_last_not_of("ACGTacgt");
			if (pos == std::string::npos) {
				m_u = Kmer(kmer);
				assert(m_g->m_bloom[m_u]);
				return;
			} else
				m_i += pos;
		}
	}

  public:
	vertex_iterator() { }

	vertex_iterator(const DBGBloom& g, size_t i)
		: m_g(&g), m_i(i)
	{
		next();
	}

	const vertex_descriptor operator*() const
	{
		return m_u;
	}

	bool operator==(const vertex_iterator& it) const
	{
		return m_i == it.m_i;
	}

	bool operator!=(const vertex_iterator& it) const
	{
		return !(*this == it);
	}

	vertex_iterator& operator++()
	{
		assert(m_i < m_g->m_fa.size());
		++m_i;
		next();
		return *this;
	}

  private:
	const DBGBloom* m_g;
	size_t m_i;
	Kmer m_u;
}; // vertex_iterator

}; // graph_traits<DBGBloom>

/** PropertyGraph vertex_index */
template<>
struct property_map<DBGBloom, vertex_index_t> {
	typedef DBGBloomIndexMap type;
	typedef type const_type;
};

} // namespace boost

// AdjacencyGraph

static inline
std::pair<graph_traits<DBGBloom>::adjacency_iterator,
	graph_traits<DBGBloom>::adjacency_iterator>
adjacent_vertices(
		graph_traits<DBGBloom>::vertex_descriptor u,
		const DBGBloom& g)
{
	typedef graph_traits<DBGBloom>::adjacency_iterator
		adjacency_iterator;
	return std::make_pair(adjacency_iterator(g, u),
			adjacency_iterator(g));
}

// IncidenceGraph

static inline
graph_traits<DBGBloom>::degree_size_type
out_degree(
		graph_traits<DBGBloom>::vertex_descriptor u,
		const DBGBloom& g)
{
	typedef graph_traits<DBGBloom>::adjacency_iterator Ait;
	std::pair<Ait, Ait> adj = adjacent_vertices(u, g);
	return std::distance(adj.first, adj.second);
}

static inline
std::pair<graph_traits<DBGBloom>::out_edge_iterator,
	graph_traits<DBGBloom>::out_edge_iterator>
out_edges(
		graph_traits<DBGBloom>::vertex_descriptor u,
		const DBGBloom& g)
{
	typedef graph_traits<DBGBloom>::out_edge_iterator Oit;
	return std::make_pair(Oit(g, u), Oit(g));
}

// BidirectionalGraph

static inline
graph_traits<DBGBloom>::degree_size_type
in_degree(graph_traits<DBGBloom>::vertex_descriptor u,
		const DBGBloom& g)
{
	return out_degree(reverseComplement(u), g);
}

static inline
std::pair<graph_traits<DBGBloom>::in_edge_iterator,
	graph_traits<DBGBloom>::in_edge_iterator>
in_edges(
		graph_traits<DBGBloom>::vertex_descriptor u,
		const DBGBloom& g)
{
	typedef graph_traits<DBGBloom>::in_edge_iterator Iit;
	return std::make_pair(Iit(g, u), Iit(g));
}

// VertexListGraph

static inline
graph_traits<DBGBloom>::vertices_size_type
num_vertices(const DBGBloom& g)
{
	return g.m_bloom.popcount();
}

static inline
std::pair<graph_traits<DBGBloom>::vertex_iterator,
	graph_traits<DBGBloom>::vertex_iterator>
vertices(const DBGBloom& g)
{
	typedef graph_traits<DBGBloom>::vertex_iterator Vit;
	return std::make_pair(Vit(g, 0), Vit(g, g.m_fa.size()));
}

// EdgeListGraph

static inline
graph_traits<DBGBloom>::edges_size_type
num_edges(const DBGBloom& g)
{
	typedef graph_traits<DBGBloom>::vertex_iterator Vit;
	size_t n = 0;
	std::pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit)
		n += out_degree(*uit, g);
	return n;
}

// PropertyGraph vertex_index

static inline
DBGBloomIndexMap
get(vertex_index_t, const DBGBloom& g)
{
	return DBGBloomIndexMap(g);
}

static inline
graph_traits<DBGBloom>::vertices_size_type
get(vertex_index_t tag, const DBGBloom& g,
		graph_traits<DBGBloom>::vertex_descriptor u)
{
	return get(get(tag, g), u);
}

// PropertyGraph

/** Return the reverse complement of the specified k-mer. */
static inline
graph_traits<DBGBloom>::vertex_descriptor
get(vertex_complement_t, const DBGBloom&,
		graph_traits<DBGBloom>::vertex_descriptor u)
{
	return reverseComplement(u);
}

/** Return the name of the specified vertex. */
static inline
Kmer get(vertex_name_t, const DBGBloom&,
		graph_traits<DBGBloom>::vertex_descriptor u)
{
	return u;
}

static inline
bool
get(vertex_removed_t, const DBGBloom&,
		graph_traits<DBGBloom>::vertex_descriptor)
{
	return false;
}

static inline
no_property
get(vertex_bundle_t, const DBGBloom&,
		graph_traits<DBGBloom>::edge_descriptor)
{
	return no_property();
}

static inline
no_property
get(edge_bundle_t, const DBGBloom&,
		graph_traits<DBGBloom>::edge_descriptor)
{
	return no_property();
}

#endif
