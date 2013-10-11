/**
 * de Bruijn Graph data structure using an FM-Index
 * Copyright 2013 Shaun Jackman
 */

#ifndef DBGFM_H
#define DBGFM_H 1

#include "Assembly/KmerData.h"
#include "Common/Kmer.h"
#include "Common/IOUtil.h"
#include "FMIndex/bit_array.h"

#include "fm_index.h"

#include <cassert>
#include <cstdlib> // for abort
#include <iostream>
#include <string>
#include <vector>

/** Mark a vertex in the out direction. */
enum vertex_mark_out_t { vertex_mark_out };

/** Mark a vertex in the in direction. */
enum vertex_mark_in_t { vertex_mark_in };

// Properties
namespace boost {
	BOOST_INSTALL_PROPERTY(vertex, mark_out);
	BOOST_INSTALL_PROPERTY(vertex, mark_in);
}

/** The index and orientation of a vertex. */
struct vertex_index_sense {
	size_t i;
	bool sense;
};

/** de Bruijn Graph data structure using an FM-Index */
class DBGFM {
  public:
	typedef KmerData vertex_bundled;

	/** The size of a k-mer */
	unsigned m_k;

	/** The FM index */
	FMIndex m_fm;

	/** Map suffix array indices to vertex indices. */
	wat_array::BitArray m_rank;

	/** Vertex property map */
	std::vector<vertex_bundled> m_vpmap;

	/** Construct a de Bruijn graph. */
	DBGFM(unsigned k, const std::string& fmPath)
		: m_k(k), m_fm(fmPath), m_rank(m_fm.getBWLen())
	{
		// Map suffix array indices to vertex indices.
		std::string s(m_k, '$');
		size_t index = 0;
		size_t n = m_fm.getBWLen() - 1;
		for (size_t i = 0; i < n; ++i) {
			char c = m_fm.getChar(index);
			assert(c != EOF);
			s.resize(s.size() - 1);
			s.insert(s.begin(), c);
			index = m_fm.LF(index);
			m_rank.SetBit(s.find('$') == std::string::npos, index);
		}
		m_rank.Build();

		m_vpmap.resize(m_rank.one_num());

#ifndef NDEBUG
		// Check that we have the expected number of k-mer.
		size_t numChars = m_fm.getBWLen() - 1;
		size_t numStrings = m_fm.getNumStrings();
		assert(numChars > numStrings * m_k);
		assert(numChars - numStrings * m_k == m_rank.one_num());
#endif
	}

	/** Return the vertex property. */
	const vertex_bundled& operator[](size_t sai) const
	{
		assert(sai < m_rank.length());
		size_t i = m_rank.Rank(1, sai);
		assert(i < m_vpmap.size());
		return m_vpmap[i];
	}

	/** Return the vertex property. */
	vertex_bundled& operator[](size_t sai)
	{
		assert(sai < m_rank.length());
		size_t i = m_rank.Rank(1, sai);
		assert(i < m_vpmap.size());
		return m_vpmap[i];
	}

	/** Return the reoriented vertex property. */
	vertex_bundled operator[](const vertex_index_sense& x) const
	{
		return x.sense ? ~(*this)[x.i] : (*this)[x.i];
	}

  private:
	/** Copy constructor. */
	DBGFM(const DBGFM&);
}; // class DBGFM

/** Return the index of this vertex and its orientation. */
static inline
vertex_index_sense
orientVertex(const DBGFM& g, const Kmer& u)
{
	struct vertex_index_sense x;
	std::pair<size_t, size_t> sai = g.m_fm.findInterval(u.str());
	if (sai.first <= sai.second) {
		x.i = sai.first;
		x.sense = false;
		return x;
	}
	sai = g.m_fm.findInterval(reverseComplement(u).str());
	if (sai.first <= sai.second) {
		x.i = sai.first;
		x.sense = true;
		return x;
	}
	std::cerr << "error: " << __func__ << ": " << u << '\n';
	assert(false);
	abort();
}

// Graph

namespace boost {

template <>
struct graph_traits<DBGFM> {
	// Graph
	typedef Kmer vertex_descriptor;
	typedef boost::directed_tag directed_category;
	struct traversal_category
		: boost::adjacency_graph_tag, boost::vertex_list_graph_tag
		{ };
	typedef boost::disallow_parallel_edge_tag edge_parallel_category;

	// IncidenceGraph
	typedef std::pair<vertex_descriptor, vertex_descriptor>
		edge_descriptor;
	typedef unsigned degree_size_type;
	typedef void out_edge_iterator;

	// BidirectionalGraph
	typedef void in_edge_iterator;

	// VertexListGraph
	typedef size_t vertices_size_type;

	// EdgeListGraph
	typedef void edge_iterator;
	typedef void edges_size_type;

// AdjacencyGraph
/** Iterate through the adjacent vertices of a vertex. */
struct adjacency_iterator
	: public std::iterator<std::input_iterator_tag, vertex_descriptor>
{
	/** Skip to the next edge that is present. */
	void next()
	{
		for (; m_i < NUM_BASES && !m_adj.checkBase(m_i); m_i++) {
		}
		if (m_i < NUM_BASES)
			m_v.setLastBase(SENSE, m_i);
	}

  public:
	adjacency_iterator() : m_i(NUM_BASES) { }

	adjacency_iterator(
			vertex_descriptor u, SeqExt adj)
		: m_v(u), m_adj(adj), m_i(0)
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
	vertex_descriptor m_v;
	SeqExt m_adj;
	short unsigned m_i;
}; // adjacency_iterator

// VertexListGraph
/** Iterate through the vertices of this graph. */
struct vertex_iterator
	: public std::iterator<std::input_iterator_tag, vertex_descriptor>
{
	typedef graph_traits<DBGFM>::vertices_size_type It;

	/** Increment */
	void increment()
	{
		assert(m_it < m_last);
		++m_it;
		char c = m_g->m_fm.getChar(m_ui);
		if (c == EOF) {
			assert(m_it == m_last);
			return;
		}
		assert(!m_s.empty());
		m_s.resize(m_s.size() - 1);
		m_s.insert(m_s.begin(), c);
		m_ui = m_g->m_fm.LF(m_ui);
	}

	/** Skip to the next vertex that is present. */
	void next()
	{
		for (; m_it < m_last; increment()) {
			if (m_s.find('$') == std::string::npos) {
				m_u = vertex_descriptor(m_s);
				break;
			}
		}
	}

  public:
	vertex_iterator(const DBGFM& g, const It& it)
		: m_g(&g), m_last(m_g->m_fm.getBWLen()), m_it(it),
		m_s(g.m_k, '$'), m_ui(0)
	{
		next();
	}

	const vertex_descriptor operator*() const
	{
		assert(m_it < m_last);
		return m_u;
	}

	bool operator==(const vertex_iterator& it) const
	{
		return m_it == it.m_it;
	}

	bool operator!=(const vertex_iterator& it) const
	{
		return !(*this == it);
	}

	vertex_iterator& operator++()
	{
		assert(m_it < m_last);
		increment();
		next();
		return *this;
	}

	/** Return the vertex index of this iterator. */
	vertices_size_type index() const
	{
		assert(m_it < m_last);
		return m_ui;
	}

	/** Return the position in the original text. */
	size_t position()
	{
		return m_last - m_it;
	}

  private:
	const DBGFM* m_g;
	It m_last;
	It m_it;
	std::string m_s;
	vertex_descriptor m_u;
	vertices_size_type m_ui;
}; // vertex_iterator

}; // graph_traits<DBGFM>

} // namespace boost

// PropertyGraph

static inline
graph_traits<DBGFM>::vertices_size_type
get(vertex_index_t, const DBGFM& g,
		graph_traits<DBGFM>::vertex_descriptor u)
{
	typedef graph_traits<DBGFM>::vertices_size_type Vi;
	std::pair<Vi, Vi> x = g.m_fm.findInterval(u.str());
	assert(x.first <= x.second);
	return x.first;
}

static inline
const vertex_bundle_type<DBGFM>::type&
get(vertex_bundle_t, const DBGFM& g,
		graph_traits<DBGFM>::vertex_descriptor u)
{
	typedef graph_traits<DBGFM>::vertices_size_type Vi;
	Vi ui = get(vertex_index, g, u);
	assert(ui < g.m_rank.length());
	return g[ui];
}

static inline
vertex_bundle_type<DBGFM>::type&
get(vertex_bundle_t, DBGFM& g,
		graph_traits<DBGFM>::vertex_descriptor u)
{
	typedef graph_traits<DBGFM>::vertices_size_type Vi;
	Vi ui = get(vertex_index, g, u);
	assert(ui < g.m_rank.length());
	return g[ui];
}

static inline
const vertex_bundle_type<DBGFM>::type&
get(vertex_bundle_t, const DBGFM& g,
		const graph_traits<DBGFM>::vertex_iterator& uit)
{
	typedef graph_traits<DBGFM>::vertices_size_type Vi;
	Vi ui = uit.index();
	assert(ui < g.m_rank.length());
	return g[ui];
}

static inline
vertex_bundle_type<DBGFM>::type&
get(vertex_bundle_t, DBGFM& g,
		const graph_traits<DBGFM>::vertex_iterator& uit)
{
	typedef graph_traits<DBGFM>::vertices_size_type Vi;
	Vi ui = uit.index();
	assert(ui < g.m_rank.length());
	return g[ui];
}

// IncidenceGraph

static inline
graph_traits<DBGFM>::degree_size_type
out_degree(
		graph_traits<DBGFM>::vertex_descriptor u,
		const DBGFM& g)
{
	return get(vertex_bundle, g, u).getExtension(SENSE).outDegree();
}

// BidirectionalGraph

static inline
graph_traits<DBGFM>::degree_size_type
in_degree(graph_traits<DBGFM>::vertex_descriptor u,
		const DBGFM& g)
{
	return get(vertex_bundle, g, u).getExtension(
			ANTISENSE).outDegree();
}

// AdjacencyGraph

static inline
std::pair<graph_traits<DBGFM>::adjacency_iterator,
	graph_traits<DBGFM>::adjacency_iterator>
adjacent_vertices(
		graph_traits<DBGFM>::vertex_descriptor u,
		const DBGFM& g)
{
	typedef graph_traits<DBGFM>::adjacency_iterator
		adjacency_iterator;
	SeqExt adj = get(vertex_bundle, g, u).getExtension(SENSE);
	return std::make_pair(adjacency_iterator(u, adj),
			adjacency_iterator());
}

// VertexListGraph

static inline
graph_traits<DBGFM>::vertices_size_type
num_vertices(const DBGFM& g)
{
	return g.m_rank.one_num();
}

static inline
std::pair<graph_traits<DBGFM>::vertex_iterator,
	graph_traits<DBGFM>::vertex_iterator>
vertices(const DBGFM& g)
{
	typedef graph_traits<DBGFM>::vertex_iterator Vit;
	return std::make_pair(Vit(g, 0), Vit(g, g.m_fm.getBWLen()));
}

// PropertyGraph

/** Return the reverse complement of the specified k-mer. */
static inline
graph_traits<DBGFM>::vertex_descriptor
get(vertex_complement_t, const DBGFM&,
		graph_traits<DBGFM>::vertex_descriptor u)
{
	return reverseComplement(u);
}

static inline
bool
get(vertex_removed_t, const DBGFM& g,
		graph_traits<DBGFM>::vertex_descriptor u)
{
	return get(vertex_bundle, g, u).deleted();
}

static inline
void
put(vertex_removed_t, DBGFM& g,
		graph_traits<DBGFM>::vertex_descriptor u, bool flag)
{
	(void)flag;
	assert(flag);
	get(vertex_bundle, g, u).setFlag(SF_DELETE);
}

static inline
void
put(vertex_mark_out_t, DBGFM& g,
		graph_traits<DBGFM>::vertex_descriptor u, bool flag)
{
	(void)flag;
	assert(flag);
	vertex_index_sense x = orientVertex(g, u);
	assert(x.i < g.m_rank.length());
	g[x.i].setFlag(
			x.sense ? SF_MARK_ANTISENSE : SF_MARK_SENSE);
}

static inline
void
put(vertex_mark_in_t, DBGFM& g,
		graph_traits<DBGFM>::vertex_descriptor u, bool flag)
{
	(void)flag;
	assert(flag);
	vertex_index_sense x = orientVertex(g, u);
	assert(x.i < g.m_rank.length());
	g[x.i].setFlag(
			x.sense ? SF_MARK_SENSE : SF_MARK_ANTISENSE);
}

static inline
no_property
get(edge_bundle_t, const DBGFM&,
		graph_traits<DBGFM>::edge_descriptor)
{
	return no_property();
}

// Subgraph

static inline
bool
vertex_exists(graph_traits<DBGFM>::vertex_descriptor u,
		const DBGFM& g)
{
	return g.m_fm.count(u.str()) > 0
		|| g.m_fm.count(get(vertex_complement, g, u).str()) > 0;
}

// Distributed graph

static inline
void
pumpNetwork(const DBGFM&)
{
}

#endif
