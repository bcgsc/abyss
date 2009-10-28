#ifndef ISEQUENCECOLLECTION_H
#define ISEQUENCECOLLECTION_H 1

#include "config.h"
#include "PackedSeq.h"

struct PackedSeqHasher
{
	size_t operator()(const PackedSeq& o) const
	{
		return o.getHashCode();
	}
};

#if HAVE_GOOGLE_SPARSE_HASH_SET
# include <google/sparse_hash_set>
typedef google::sparse_hash_set<PackedSeq,
		PackedSeqHasher> SequenceDataHash;
#else
# include "HashSet.h"
typedef hash_set<PackedSeq, PackedSeqHasher> SequenceDataHash;
#endif

// Interface class for a sequence collection (the lowest level of storage of a large number of sequences)
// This pure virtual class defines the minimum set of functions a sequence collection must provide
class ISequenceCollection
{
	public:
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

		// Return the number of sequences in the collection
		virtual size_t count() const = 0;

		virtual void printLoad() const = 0;

		// Clear the specified flag from every sequence in the collection
		virtual void wipeFlag(SeqFlag flag) = 0;

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

		virtual bool getSeqData(const Kmer& seq, ExtensionRecord& ext,
				int& multiplicity) const = 0;
		virtual const PackedSeq& getSeqAndData(
				const Kmer& key) const = 0;

		// Receive and dispatch packets if necessary.
		virtual unsigned pumpNetwork() = 0;

		virtual iterator begin() = 0;
		virtual const_iterator begin() const = 0;
		virtual iterator end() = 0;
		virtual const_iterator end() const = 0;

		// Observer pattern
		typedef void (*SeqObserver)(ISequenceCollection* c,
				const PackedSeq& seq);
		virtual void attach(SeqObserver f) = 0;
		virtual void detach(SeqObserver f) = 0;

		virtual void load(const char *path) = 0;

		virtual void setColourSpace(bool flag) = 0;
};

#endif
