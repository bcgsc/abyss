#ifndef PAIREDDBG_ISEQUENCECOLLECTION_H
#define PAIREDDBG_ISEQUENCECOLLECTION_H 1

#include "config.h"
#include "KmerPair.h"
#include "KmerPairData.h"

#if HAVE_GOOGLE_SPARSE_HASH_MAP
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<KmerPair, KmerPairData, hash<KmerPair> >
	SequenceDataHash;
#else
# include "UnorderedMap.h"
typedef unordered_map<KmerPair, KmerPairData, hash<KmerPair> >
	SequenceDataHash;
#endif

/** A hash table mapping vertices to vertex properties. */
class ISequenceCollection
{
	public:
		typedef SequenceDataHash::key_type key_type;
		typedef SequenceDataHash::mapped_type mapped_type;
		typedef SequenceDataHash::value_type value_type;
		typedef SequenceDataHash::iterator iterator;
		typedef SequenceDataHash::const_iterator const_iterator;
		typedef Dinuc OutEdgeDescriptor;
		typedef SeqExt EdgeSet;

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
				extDirection dir, EdgeSet ext) = 0;

		/** Remove the specified edge of this k-mer. */
		void removeExtension(const key_type& seq,
				extDirection dir, OutEdgeDescriptor base)
		{
			removeExtension(seq, dir, EdgeSet(base));
		}

		/** Remove all the edges of this k-mer. */
		void clearExtensions(const key_type& seq, extDirection dir)
		{
			removeExtension(seq, dir, EdgeSet::mask(0xf));
		}

		virtual bool setBaseExtension(const key_type& seq,
				extDirection dir, OutEdgeDescriptor base) = 0;

		// Receive and dispatch packets if necessary.
		virtual size_t pumpNetwork() = 0;

		virtual iterator begin() = 0;
		virtual const_iterator begin() const = 0;
		virtual iterator end() = 0;
		virtual const_iterator end() const = 0;

		// Observer pattern
		typedef void (*SeqObserver)(ISequenceCollection* c,
				const value_type& seq);
		virtual void attach(SeqObserver f) = 0;
		virtual void detach(SeqObserver f) = 0;

		virtual void load(const char *path) = 0;

		virtual void setColourSpace(bool flag) = 0;
};

#endif
