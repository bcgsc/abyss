#ifndef ASSEMBLY_DBG_H
#define ASSEMBLY_DBG_H 1

#include "config.h"
#include "Assembly/Options.h"
#include "Common/Log.h"
#include "Common/MemoryUtil.h"
#include "Common/Options.h"
#include "Common/StringUtil.h" // for toSI
#include "Common/Timer.h"
#include "Graph/Properties.h"

#include <boost/graph/graph_traits.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <utility>

using boost::graph_traits;

/** A hash table mapping vertices to vertex properties. */
class SequenceCollectionHash
{
	public:
		typedef SequenceDataHash::key_type key_type;
		typedef SequenceDataHash::mapped_type mapped_type;
		typedef SequenceDataHash::value_type value_type;
		typedef SequenceDataHash::iterator iterator;
		typedef SequenceDataHash::const_iterator const_iterator;

		typedef mapped_type::Symbol Symbol;
		typedef mapped_type::SymbolSet SymbolSet;
		typedef mapped_type::SymbolSetPair SymbolSetPair;

		typedef key_type vertex_descriptor;
		typedef mapped_type vertex_bundled;
		typedef std::pair<key_type, key_type> edge_descriptor;

		/** Remove the specified sequence if it exists. */
		void remove(const key_type& seq)
		{
			setFlag(seq, SF_DELETE);
		}

		/** Shrink the hash table. */
		void shrink() {
			m_data.rehash(0);
			printLoad();
		}

		/** Return the data associated with the specified key. */
		const mapped_type operator[](const key_type& key) const
		{
			bool rc;
			const_iterator it = find(key, rc);
			assert(it != m_data.end());
			return rc ? ~it->second : it->second;
		}

		iterator begin() { return m_data.begin(); }
		const_iterator begin() const { return m_data.begin(); }
		iterator end() { return m_data.end(); }
		const_iterator end() const { return m_data.end(); }

		/** Return true if this collection is empty. */
		bool empty() const { return m_data.empty(); }

		/** Return the number of sequences in this collection. */
		size_t size() const { return m_data.size(); }

		// Not a network sequence collection. Nothing to do.
		size_t pumpNetwork() { return 0; }

		/** The observer callback function. */
		typedef void (*SeqObserver)(SequenceCollectionHash* c,
				const value_type& seq);

		/** Attach the specified observer. */
		void attach(SeqObserver f)
		{
			assert(m_seqObserver == NULL);
			m_seqObserver = f;
		}

		/** Detach the specified observer. */
		void detach(SeqObserver f)
		{
			assert(m_seqObserver == f);
			(void)f;
			m_seqObserver = NULL;
		}

		bool isAdjacencyLoaded() const { return m_adjacencyLoaded; }

SequenceCollectionHash()
	: m_seqObserver(NULL), m_adjacencyLoaded(false)
{
#if HAVE_GOOGLE_SPARSE_HASH_MAP
	// sparse_hash_set uses 2.67 bits per element on a 64-bit
	// architecture and 2 bits per element on a 32-bit architecture.
	// The number of elements is rounded up to a power of two.
	if (opt::rank >= 0) {
		// Initialize sparsehash size to 2^30 (~1 billion) empty buckets.
		// Setting the initial sparsehash size to a large
		// number avoids the prohibitively slow step of resizing
		// the hash table and rehashing *every* element when the
		// maximum load factor is exceeded.
		//
		// The initial memory footprint per CPU core is
		// 2^30 * 2.67 / 8 ~= 0.358 GB. The value of 2^30 buckets
		// was chosen to accomodate a k=144 32-thread assembly of human
		// (NA24143) uncorrected Illumina reads with ~70X coverage and
		// 20,317,980,431 distinct 144-mers, without requiring sparsehash
		// resizing. For further details on the test dataset, see:
		// "ABySS 2.0: Resource-efficient assembly of large genomes
		// using a Bloom filter".
		m_data.rehash((size_t)pow(2, 30));
		m_data.min_load_factor(0.2);
	} else {
		// Allocate a big hash for a single processor.
		m_data.rehash(1<<29);
		m_data.max_load_factor(0.4);
	}
#endif
}

/** sparse_hash_set requires that set_deleted_key()
 * is called before calling erase(). This key cannot
 * be an existing kmer in m_data. This function sets
 * the deleted key and should be called after all
 * data has been loaded.
 */
void setDeletedKey()
{
#if HAVE_GOOGLE_SPARSE_HASH_MAP
	for (SequenceDataHash::iterator it = m_data.begin();
			it != m_data.end(); it++) {
		key_type rc(reverseComplement(it->first));
		bool isrc;
		SequenceDataHash::iterator search = find(rc, isrc);
		// If this is false, we should have a palindrome or we're
		// doing a SS assembly.
		if (isrc || search == m_data.end()) {
			m_data.set_deleted_key(rc);
			return;
		}
	}
	logger(1) << "error: unable to set deleted key.\n";
	exit(EXIT_FAILURE);
#else
	return;
#endif
}

/** Add the specified k-mer to this collection. */
void add(const key_type& seq, unsigned coverage = 1)
{
	bool rc;
	iterator it = find(seq, rc);
	if (it == m_data.end()) {
		m_data.insert(std::make_pair(seq, mapped_type(SENSE, coverage)));
	} else if (coverage > 0) {
		assert(!rc || !opt::ss);
		it->second.addMultiplicity(rc ? ANTISENSE : SENSE, coverage);
	}
}

/** Clean up by erasing sequences flagged as deleted.
 * @return the number of sequences erased
 */
size_t cleanup()
{
	Timer(__func__);
	size_t count = 0;
	for (iterator it = m_data.begin(); it != m_data.end();) {
		if (it->second.deleted()) {
			m_data.erase(it++);
			count++;
		} else
			++it;
	}
	shrink();
	return count;
}

/** Add an edge to this k-mer. */
bool setBaseExtension(
		const key_type& kmer, extDirection dir, Symbol base)
{
	bool rc;
	iterator it = find(kmer, rc);
	if (it == m_data.end())
		return false;
	if (opt::ss) {
		assert(!rc);
		it->second.setBaseExtension(dir, base);
	} else {
		bool palindrome = kmer.isPalindrome();
		if (!rc || palindrome)
			it->second.setBaseExtension(dir, base);
		if (rc || palindrome)
			it->second.setBaseExtension(!dir, reverseComplement(base));
	}
	return true;
}

/** Remove the specified extensions from this k-mer. */
void removeExtension(const key_type& kmer,
		extDirection dir, SymbolSet ext)
{
	bool rc;
	iterator it = find(kmer, rc);
	assert(it != m_data.end());
	if (opt::ss) {
		assert(!rc);
		it->second.removeExtension(dir, ext);
	} else {
		bool palindrome = kmer.isPalindrome();
		if (!rc || palindrome)
			it->second.removeExtension(dir, ext);
		if (rc || palindrome)
			it->second.removeExtension(!dir, ext.complement());
	}
	notify(*it);
}

/** Remove the specified edge of this vertex. */
void removeExtension(const key_type& seq, extDirection dir, Symbol base)
{
	removeExtension(seq, dir, SymbolSet(base));
}

void setFlag(const key_type& key, SeqFlag flag)
{
	bool rc;
	iterator it = find(key, rc);
	assert(it != m_data.end());
	it->second.setFlag(rc ? complement(flag) : flag);
}

/** Mark the specified sequence in both directions. */
void mark(const key_type& seq)
{
	setFlag(seq, SeqFlag(SF_MARK_SENSE | SF_MARK_ANTISENSE));
}

/** Mark the specified sequence. */
void mark(const key_type& seq, extDirection sense)
{
	setFlag(seq, sense == SENSE
			? SF_MARK_SENSE : SF_MARK_ANTISENSE);
}

/** Clear the specified flag for all vertices. */
void wipeFlag(SeqFlag flag)
{
	for (iterator it = m_data.begin();
			it != m_data.end(); ++it)
		it->second.clearFlag(flag);
}

/** Print the load of the hash table. */
void printLoad() const
{
	size_t size = m_data.size();
	size_t buckets = m_data.bucket_count();
	logger(1) << "Hash load: " << size << " / " << buckets << " = "
		<< std::setprecision(3) << (float)size / buckets
		<< " using " << toSI(getMemoryUsage()) << "B" << std::endl;
}

private:

const_iterator
find(const key_type& key) const
{
	return m_data.find(key);
}

iterator
find(const key_type& key)
{
	return m_data.find(key);
}

/** Return an iterator pointing to the specified k-mer or its
 * reverse complement. Return in rc whether the sequence is reversed.
 */
iterator
find(const key_type& key, bool& rc)
{
	iterator it = find(key);
	if (opt::ss || it != m_data.end()) {
		rc = false;
		return it;
	} else {
		rc = true;
		return find(reverseComplement(key));
	}
}

public:

/** Return an iterator pointing to the specified k-mer or its
 * reverse complement. Return in rc whether the sequence is reversed.
 */
const_iterator
find(const key_type& key, bool& rc) const
{
	const_iterator it = find(key);
	if (opt::ss || it != m_data.end()) {
		rc = false;
		return it;
	} else {
		rc = true;
		return find(reverseComplement(key));
	}
}

/** Return the sequence and data of the specified key.
 * The key sequence may not contain data. The returned sequence will
 * contain data.
 */
const value_type&
getSeqAndData(const key_type& key) const
{
	bool rc;
	const_iterator it = find(key, rc);
	// rc should not be ignored. This seems quite dubious.
	// The edges of this k-mer should be complemented.
	assert(it != m_data.end());
	return *it;
}

/** Return the data of the specified key. */
bool getSeqData(const key_type& key,
		SymbolSetPair& extRecord, int& multiplicity) const
{
	bool rc;
	const_iterator it = find(key, rc);
	assert(!rc || !opt::ss);
	if (it == m_data.end())
		return false;
	const mapped_type data = it->second;
	extRecord = rc ? data.extension().complement() : data.extension();
	multiplicity = data.getMultiplicity();
	return true;
}

/** Write this collection to disk.
 * @param path does not include the extension
 */
void store(const char* path)
{
	assert(path != NULL);
#if HAVE_GOOGLE_SPARSE_HASH_MAP
	std::ostringstream s;
	s << path;
	if (opt::rank >= 0)
		s << '-' << std::setfill('0') << std::setw(3) << opt::rank;
	s << ".kmer";
	FILE* f = fopen(s.str().c_str(), "w");
	if (f == NULL) {
		perror(s.str().c_str());
		exit(EXIT_FAILURE);
	}
	shrink();
	m_data.write_metadata(f);
	m_data.write_nopointer_data(f);
	fclose(f);
#else
	// Not supported.
	assert(false);
	exit(EXIT_FAILURE);
#endif
}

/** Load this collection from disk. */
void load(const char* path)
{
#if HAVE_GOOGLE_SPARSE_HASH_MAP
	FILE* f = fopen(path, "r");
	if (f == NULL) {
		perror(path);
		exit(EXIT_FAILURE);
	}
	m_data.read_metadata(f);
	m_data.read_nopointer_data(f);
	fclose(f);
	m_adjacencyLoaded = true;
#else
	(void)path;
	// Not supported.
	assert(false);
	exit(EXIT_FAILURE);
#endif
}

/** Indicate that this is a colour-space collection. */
void setColourSpace(bool flag)
{
	if (!m_data.empty())
		assert(opt::colourSpace == flag);
	opt::colourSpace = flag;
}

	private:
		/** Call the observers of the specified sequence. */
		void notify(const value_type& seq)
		{
			if (m_seqObserver != NULL)
				m_seqObserver(this, seq);
		}

		/** The underlying collection. */
		SequenceDataHash m_data;

		/** The observers. Only a single observer is implemented.*/
		SeqObserver m_seqObserver;

		/** Whether adjacency information has been loaded. */
		bool m_adjacencyLoaded;
};

// Forward declaration
class DBGEdgeIterator;

// Graph

namespace boost {

template <>
struct graph_traits<SequenceCollectionHash> {
	// Graph
	typedef SequenceCollectionHash Graph;
	typedef Graph::key_type vertex_descriptor;
	typedef boost::directed_tag directed_category;
	struct traversal_category
		: boost::adjacency_graph_tag, boost::vertex_list_graph_tag
		{ };
	typedef boost::disallow_parallel_edge_tag edge_parallel_category;

	// IncidenceGraph
	typedef std::pair<vertex_descriptor, vertex_descriptor>
		edge_descriptor;
	typedef unsigned degree_size_type;

	// VertexListGraph
	typedef size_t vertices_size_type;

	// EdgeListGraph
	typedef size_t edges_size_type;
	typedef DBGEdgeIterator edge_iterator;

	// Other
	typedef Graph::Symbol Symbol;
	typedef Graph::SymbolSet SymbolSet;
	static const unsigned NUM_SYMBOLS = SymbolSet::NUM;

// AdjacencyGraph
/** Iterate through the adjacent vertices of a vertex. */
struct adjacency_iterator
	: public std::iterator<std::input_iterator_tag, vertex_descriptor>
{
	/** Skip to the next edge that is present. */
	void next()
	{
		for (; m_i < NUM_SYMBOLS && !m_adj.checkBase(Symbol(m_i)); ++m_i) {
		}
		if (m_i < NUM_SYMBOLS)
			m_v.setLastBase(SENSE, Symbol(m_i));
	}

  public:
	adjacency_iterator() : m_i(NUM_SYMBOLS) { }

	adjacency_iterator(
			vertex_descriptor u, SymbolSet adj)
		: m_v(u), m_adj(adj), m_i(0)
	{
		m_v.shift(SENSE);
		next();
	}

	const vertex_descriptor& operator*() const
	{
		assert(m_i < NUM_SYMBOLS);
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
		assert(m_i < NUM_SYMBOLS);
		++m_i;
		next();
		return *this;
	}

  private:
	vertex_descriptor m_v;
	SymbolSet m_adj;
	short unsigned m_i;
}; // adjacency_iterator

// IncidenceGraph

/** Iterate through the out edges of a vertex. */
struct out_edge_iterator
	: public std::iterator<std::input_iterator_tag, edge_descriptor>
{
	/** Skip to the next edge that is present. */
	void next()
	{
		for (; m_i < NUM_SYMBOLS && !m_adj.checkBase(Symbol(m_i)); ++m_i) {
		}
		if (m_i < NUM_SYMBOLS)
			m_e.second.setLastBase(SENSE, Symbol(m_i));
	}

  public:
	out_edge_iterator() : m_i(NUM_SYMBOLS) { }

	out_edge_iterator(
			vertex_descriptor u, SymbolSet adj)
		: m_e(u, u), m_adj(adj), m_i(0)
	{
		m_e.second.shift(SENSE);
		next();
	}

	const edge_descriptor& operator*() const
	{
		assert(m_i < NUM_SYMBOLS);
		return m_e;
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
		assert(m_i < NUM_SYMBOLS);
		++m_i;
		next();
		return *this;
	}

  private:
	edge_descriptor m_e;
	SymbolSet m_adj;
	short unsigned m_i;
}; // out_edge_iterator

// BidirectionalGraph

/** Iterate through the in-edges of a vertex. */
struct in_edge_iterator
	: public std::iterator<std::input_iterator_tag, edge_descriptor>
{
	/** Skip to the next edge that is present. */
	void next()
	{
		for (; m_i < NUM_SYMBOLS && !m_adj.checkBase(Symbol(m_i)); ++m_i) {
		}
		if (m_i < NUM_SYMBOLS)
			m_e.first.setLastBase(ANTISENSE, Symbol(m_i));
	}

  public:
	in_edge_iterator() : m_i(NUM_SYMBOLS) { }

	in_edge_iterator(
			vertex_descriptor u, SymbolSet adj)
		: m_e(u, u), m_adj(adj), m_i(0)
	{
		m_e.first.shift(ANTISENSE);
		next();
	}

	const edge_descriptor& operator*() const
	{
		assert(m_i < NUM_SYMBOLS);
		return m_e;
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
		assert(m_i < NUM_SYMBOLS);
		++m_i;
		next();
		return *this;
	}

  private:
	edge_descriptor m_e;
	SymbolSet m_adj;
	short unsigned m_i;
}; // in_edge_iterator

// VertexListGraph
/** Iterate through the vertices of this graph. */
struct vertex_iterator
	: public std::iterator<std::input_iterator_tag, vertex_descriptor>
{
	typedef Graph::const_iterator It;

  public:
	vertex_iterator(const It& it) : m_it(it), m_sense(false) { }

	const vertex_descriptor operator*() const
	{
		return m_sense ? reverseComplement(m_it->first) : m_it->first;
	}

	bool operator==(const vertex_iterator& it) const
	{
		return m_it == it.m_it && m_sense == it.m_sense;
	}

	bool operator!=(const vertex_iterator& it) const
	{
		return !(*this == it);
	}

	vertex_iterator& operator++()
	{
		if (m_sense) {
			++m_it;
			m_sense = false;
		} else
			m_sense = true;
		return *this;
	}

  private:
	It m_it;
	bool m_sense;
}; // vertex_iterator

}; // graph_traits<SequenceCollectionHash>

} // namespace boost

// IncidenceGraph

static inline
std::pair<
	graph_traits<SequenceCollectionHash>::out_edge_iterator,
	graph_traits<SequenceCollectionHash>::out_edge_iterator>
out_edges(
		graph_traits<SequenceCollectionHash>::vertex_descriptor u,
		const SequenceCollectionHash& g)
{
	typedef graph_traits<SequenceCollectionHash> GTraits;
	typedef GTraits::out_edge_iterator out_edge_iterator;
	typedef GTraits::SymbolSet SymbolSet;
	SymbolSet adj = g[u].getExtension(SENSE);
	return std::make_pair(
			out_edge_iterator(u, adj),
			out_edge_iterator());
}

static inline
graph_traits<SequenceCollectionHash>::degree_size_type
out_degree(
		graph_traits<SequenceCollectionHash>::vertex_descriptor u,
		const SequenceCollectionHash& g)
{
	return g[u].getExtension(SENSE).outDegree();
}

// BidirectionalGraph

static inline
std::pair<
	graph_traits<SequenceCollectionHash>::in_edge_iterator,
	graph_traits<SequenceCollectionHash>::in_edge_iterator>
in_edges(
		graph_traits<SequenceCollectionHash>::vertex_descriptor u,
		const SequenceCollectionHash& g)
{
	typedef graph_traits<SequenceCollectionHash> GTraits;
	typedef GTraits::in_edge_iterator in_edge_iterator;
	typedef GTraits::SymbolSet SymbolSet;
	SymbolSet adj = g[u].getExtension(ANTISENSE);
	return std::make_pair(
			in_edge_iterator(u, adj),
			in_edge_iterator());
}

static inline
graph_traits<SequenceCollectionHash>::degree_size_type
in_degree(graph_traits<SequenceCollectionHash>::vertex_descriptor u,
		const SequenceCollectionHash& g)
{
	return g[u].getExtension(ANTISENSE).outDegree();
}

// AdjacencyGraph

static inline
std::pair<graph_traits<SequenceCollectionHash>::adjacency_iterator,
	graph_traits<SequenceCollectionHash>::adjacency_iterator>
adjacent_vertices(
		graph_traits<SequenceCollectionHash>::vertex_descriptor u,
		const SequenceCollectionHash& g)
{
	typedef graph_traits<SequenceCollectionHash>::adjacency_iterator
		adjacency_iterator;
	typedef graph_traits<SequenceCollectionHash>::SymbolSet SymbolSet;
	SymbolSet adj = g[u].getExtension(SENSE);
	return std::make_pair(adjacency_iterator(u, adj),
			adjacency_iterator());
}

// VertexListGraph

static inline
std::pair<graph_traits<SequenceCollectionHash>::vertex_iterator,
	graph_traits<SequenceCollectionHash>::vertex_iterator>
vertices(const SequenceCollectionHash& g)
{
	return std::make_pair(g.begin(), g.end());
}

// EdgeListGraph

/** Iterate through the edges of this graph. */
class DBGEdgeIterator
	: public std::iterator<std::input_iterator_tag,
		graph_traits<SequenceCollectionHash>::edge_descriptor>
{
	typedef graph_traits<SequenceCollectionHash> GTraits;
	typedef GTraits::adjacency_iterator adjacency_iterator;
	typedef GTraits::edge_descriptor edge_descriptor;
	typedef GTraits::edge_iterator edge_iterator;
	typedef GTraits::vertex_iterator vertex_iterator;

	void nextVertex()
	{
		vertex_iterator vlast = vertices(*m_g).second;
		for (; m_vit != vlast; ++m_vit) {
			std::pair<adjacency_iterator, adjacency_iterator>
				adj = adjacent_vertices(*m_vit, *m_g);
			if (adj.first != adj.second) {
				m_eit = adj.first;
				return;
			}
		}
		// Set m_eit to a known value.
		static const adjacency_iterator s_eitNULL;
		m_eit = s_eitNULL;
	}

  public:
	DBGEdgeIterator(const SequenceCollectionHash* g, const vertex_iterator& vit)
		: m_g(g), m_vit(vit)
	{
		nextVertex();
	}

	edge_descriptor operator*() const
	{
		return edge_descriptor(*m_vit, *m_eit);
	}

	bool operator==(const edge_iterator& it) const
	{
		return m_vit == it.m_vit && m_eit == it.m_eit;
	}

	bool operator!=(const edge_iterator& it) const
	{
		return !(*this == it);
	}

	edge_iterator& operator++()
	{
		if (++m_eit == adjacent_vertices(*m_vit, *m_g).second) {
			++m_vit;
			nextVertex();
		}
		return *this;
	}

	edge_iterator operator++(int)
	{
		edge_iterator it = *this;
		++*this;
		return it;
	}

  private:
	const SequenceCollectionHash* m_g;
	vertex_iterator m_vit;
	adjacency_iterator m_eit;
}; // DBGEdgeIterator

/** Iterate through the edges of this graph. */
static inline
std::pair<
	graph_traits<SequenceCollectionHash>::edge_iterator,
	graph_traits<SequenceCollectionHash>::edge_iterator>
edges(const SequenceCollectionHash& g)
{
	
	typedef graph_traits<SequenceCollectionHash> GTraits;
	typedef GTraits::vertex_iterator vertex_iterator;
	typedef GTraits::edge_iterator edge_iterator;
	std::pair<vertex_iterator, vertex_iterator> uit = vertices(g);
	return std::make_pair(
			edge_iterator(&g, uit.first),
			edge_iterator(&g, uit.second));
}

// EdgeMutableGraph

/** Remove the edge (u,v) from the graph. */
static inline
void
remove_edge(
	graph_traits<SequenceCollectionHash>::vertex_descriptor u,
	graph_traits<SequenceCollectionHash>::vertex_descriptor v,
	SequenceCollectionHash& g)
{
	g.removeExtension(u, SENSE, v.back());
}

/** Remove the edge e from the graph. */
static inline
void
remove_edge(
		graph_traits<SequenceCollectionHash>::edge_descriptor e,
		SequenceCollectionHash& g)
{
	remove_edge(source(e, g), target(e, g), g);
}

// PropertyGraph

/** Return the reverse complement of the specified k-mer. */
static inline
graph_traits<SequenceCollectionHash>::vertex_descriptor
get(vertex_complement_t, const SequenceCollectionHash&,
		graph_traits<SequenceCollectionHash>::vertex_descriptor u)
{
	return reverseComplement(u);
}

/** Return whether this vertex has been removed. */
static inline
bool get(vertex_removed_t, const SequenceCollectionHash& g,
		graph_traits<SequenceCollectionHash>::vertex_descriptor u)
{
	return g.getSeqAndData(u).second.deleted();
}

/** Return the name of this vertex. */
static inline
std::string
get(vertex_name_t, const SequenceCollectionHash&,
		graph_traits<SequenceCollectionHash>::vertex_descriptor u)
{
	return u.str();
}

/** Return the properties of this vertex. */
static inline
vertex_bundle_type<SequenceCollectionHash>::type
get(vertex_bundle_t, const SequenceCollectionHash& g,
		graph_traits<SequenceCollectionHash>::vertex_descriptor u)
{
	return g[u];
}

/** Return the coverage of this vertex. */
static inline
unsigned
get(vertex_coverage_t, const SequenceCollectionHash& g,
		graph_traits<SequenceCollectionHash>::vertex_descriptor u)
{
	return g[u].getMultiplicity();
}

/** Return the properties of this edge. */
static inline
no_property get(edge_bundle_t, const SequenceCollectionHash&,
		graph_traits<SequenceCollectionHash>::edge_descriptor)
{
	return no_property();
}

#endif
