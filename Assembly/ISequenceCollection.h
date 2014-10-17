#ifndef ASSEMBLY_ISEQUENCECOLLECTION_H
#define ASSEMBLY_ISEQUENCECOLLECTION_H 1

#include <boost/graph/graph_traits.hpp>

/** A hash table mapping vertices to vertex properties. */
class ISequenceCollection
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

		typedef boost::directed_tag directed_category;
		typedef boost::disallow_parallel_edge_tag edge_parallel_category;
		struct traversal_category
			: boost::adjacency_graph_tag, boost::vertex_list_graph_tag
			{ };

		virtual ~ISequenceCollection() { }

		virtual void add(const key_type& seq, unsigned coverage = 1) = 0;
		virtual void remove(const key_type& seq) = 0;

		virtual void setFlag(const key_type& seq, SeqFlag flag) = 0;

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

		virtual bool empty() const = 0;

		virtual void printLoad() const = 0;

		virtual void removeExtension(const key_type& seq,
				extDirection dir, SymbolSet ext) = 0;

		virtual bool setBaseExtension(const key_type& seq,
				extDirection dir, Symbol base) = 0;

		// Receive and dispatch packets if necessary.
		virtual size_t pumpNetwork() = 0;

		virtual iterator begin() = 0;
		virtual const_iterator begin() const = 0;
		virtual iterator end() = 0;
		virtual const_iterator end() const = 0;

		virtual void load(const char *path) = 0;

		virtual void setColourSpace(bool flag) = 0;
};

#endif
