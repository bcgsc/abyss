#ifndef ISEQUENCECOLLECTION_H
#define ISEQUENCECOLLECTION_H 1

#include "config.h"
#include "Kmer.h"
#include "KmerData.h"

#if HAVE_GOOGLE_SPARSE_HASH_MAP
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<Kmer, KmerData,
		hashKmer> SequenceDataHash;
#else
# include "HashMap.h"
typedef hash_map<Kmer, KmerData, hashKmer> SequenceDataHash;
#endif

/** The interface of a map of Kmer to KmerData. */
class ISequenceCollection
{
	public:
		typedef SequenceDataHash::value_type value_type;
		typedef SequenceDataHash::iterator iterator;
		typedef SequenceDataHash::const_iterator const_iterator;

		virtual ~ISequenceCollection() { }

		virtual void add(const Kmer& seq) = 0;
		virtual void remove(const Kmer& seq) = 0;

		virtual void setFlag(const Kmer& seq, SeqFlag flag) = 0;

		/** Mark the specified sequence. */
		void mark(const Kmer& seq, extDirection sense = SENSE)
		{
			setFlag(seq, sense == SENSE
					? SF_MARK_SENSE : SF_MARK_ANTISENSE);
		}

		virtual bool empty() const = 0;

		virtual void printLoad() const = 0;

		virtual void removeExtension(const Kmer& seq,
				extDirection dir, SeqExt ext) = 0;

		/** Remove the specified edge of this k-mer. */
		void removeExtension(const Kmer& seq,
				extDirection dir, uint8_t base)
		{
			removeExtension(seq, dir, SeqExt(base));
		}

		/** Remove all the edges of this k-mer. */
		void clearExtensions(const Kmer& seq, extDirection dir)
		{
			removeExtension(seq, dir, SeqExt::mask(0xf));
		}

		virtual bool setBaseExtension(const Kmer& seq,
				extDirection dir, uint8_t base) = 0;

		// Receive and dispatch packets if necessary.
		virtual unsigned pumpNetwork() = 0;

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
