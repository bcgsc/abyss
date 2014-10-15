#ifndef SEQUENCECOLLECTION_H
#define SEQUENCECOLLECTION_H 1

#include "Graph/Properties.h"
#include "ISequenceCollection.h"
#include <boost/graph/graph_traits.hpp>
#include <cassert>
#include <utility>

using boost::graph_traits;

/** A hash table mapping vertices to vertex properties. */
class SequenceCollectionHash : public ISequenceCollection
{
	public:
		typedef SequenceDataHash::key_type key_type;
		typedef SequenceDataHash::mapped_type mapped_type;
		typedef SequenceDataHash::value_type value_type;
		typedef mapped_type vertex_property_type;
		typedef mapped_type vertex_bundled;
		typedef no_property edge_property_type;
		typedef no_property edge_bundled;

		SequenceCollectionHash();

		void add(const key_type& seq, unsigned coverage = 1);

		/** Remove the specified sequence if it exists. */
		void remove(const key_type& seq)
		{
			setFlag(seq, SF_DELETE);
		}

		// Clean up by erasing sequences flagged as deleted.
		size_t cleanup();

		/** Shrink the hash table. */
		void shrink() {
			m_data.rehash(0);
			printLoad();
		}

		// Print the load of the hash table.
		void printLoad() const;

		// Set flag for sequence seq
		void setFlag(const key_type& seq, SeqFlag flag);

		// Clear the specified flag from every sequence in the
		// collection.
		void wipeFlag(SeqFlag flag);

		bool setBaseExtension(const key_type& seq, extDirection dir,
				uint8_t base);
		void removeExtension(const key_type& seq,
				extDirection dir, SeqExt ext);

		// get the extensions of a sequence
		bool getSeqData(const key_type& seq,
				ExtensionRecord& extRecord, int& multiplicity) const;

		const value_type& getSeqAndData(const key_type& key) const;

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

		void load(const char *path);
		void store(const char* path);
		bool isAdjacencyLoaded() const { return m_adjacencyLoaded; }
		void setColourSpace(bool flag);
		void setDeletedKey();

	private:
		iterator find(const key_type& key) { return m_data.find(key); }
		const_iterator find(const key_type& key) const
		{
			return m_data.find(key);
		}

		iterator find(const key_type& key, bool& rc);
		const_iterator find(const key_type& key, bool& rc) const;

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

// Graph

namespace boost {

template <>
struct graph_traits<SequenceCollectionHash> {
	// Graph
	typedef SequenceCollectionHash::key_type vertex_descriptor;
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
	typedef SequenceCollectionHash::const_iterator It;

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
graph_traits<SequenceCollectionHash>::degree_size_type
out_degree(
		graph_traits<SequenceCollectionHash>::vertex_descriptor u,
		const SequenceCollectionHash& g)
{
	return g[u].getExtension(SENSE).outDegree();
}

// BidirectionalGraph

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
	SeqExt adj = g[u].getExtension(SENSE);
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

// PropertyGraph

/** Return the reverse complement of the specified k-mer. */
static inline
graph_traits<SequenceCollectionHash>::vertex_descriptor
get(vertex_complement_t, const SequenceCollectionHash&,
		graph_traits<SequenceCollectionHash>::vertex_descriptor u)
{
	return reverseComplement(u);
}

static inline
bool get(vertex_removed_t, const SequenceCollectionHash& g,
		graph_traits<SequenceCollectionHash>::vertex_descriptor u)
{
	return g.getSeqAndData(u).second.deleted();
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

static inline
no_property get(edge_bundle_t, const SequenceCollectionHash&,
		graph_traits<SequenceCollectionHash>::edge_descriptor)
{
	return no_property();
}

#endif
