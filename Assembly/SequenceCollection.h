#ifndef SEQUENCECOLLECTION_H
#define SEQUENCECOLLECTION_H 1

#include "ISequenceCollection.h"
#include "PackedSeq.h"
#include <cassert>

/** A k-mer graph. A map of k-mer to edges. */
class SequenceCollectionHash : public ISequenceCollection
{
	public:
		SequenceCollectionHash();
		~SequenceCollectionHash();

		void add(const Kmer& seq);

		/** Remove the specified sequence if it exists. */
		void remove(const Kmer& seq)
		{
			setFlag(seq, SF_DELETE);
		}

		// Clean up by erasing sequences flagged as deleted.
		unsigned cleanup();

		// Print the load of the hash table.
		void printLoad() const;

		// Set flag for sequence seq
		void setFlag(const Kmer& seq, SeqFlag flag);

		// Clear the specified flag from every sequence in the collection
		void wipeFlag(SeqFlag flag);

		bool setBaseExtension(const Kmer& seq, extDirection dir,
				uint8_t base);
		void removeExtension(const Kmer& seq,
				extDirection dir, SeqExt ext);

		// get the extensions of a sequence
		bool getSeqData(const Kmer& seq,
				ExtensionRecord& extRecord, int& multiplicity) const;

		const PackedSeq& getSeqAndData(const Kmer& key) const;

		iterator begin() { return m_pSequences->begin(); }
		const_iterator begin() const { return m_pSequences->begin(); }
		iterator end() { return m_pSequences->end(); }
		const_iterator end() const { return m_pSequences->end(); }

		/** Return the number of sequences in this collection. */
		size_t count() const { return m_pSequences->size(); }

		// Not a network sequence collection. Nothing to do.
		unsigned pumpNetwork() { return 0; }

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
			m_seqObserver = NULL;
		}

		void load(const char *path);
		void store(const char* path) const;
		bool isAdjacencyLoaded() const { return m_adjacencyLoaded; }
		void setColourSpace(bool flag);

	private:
		typedef std::pair<SequenceDataHash::iterator,
				SequenceDataHash::iterator> iteratorPair;

		iteratorPair findBoth(const Kmer& seq);
		iterator find(const Kmer& key);
		const_iterator find(const Kmer& key) const;
		iterator find(const Kmer& key, bool& rc);
		const_iterator find(const Kmer& key, bool& rc) const;
		const PackedSeq& getSeqAndData(
				const iteratorPair& iters) const;

		bool setBaseExtensionByIter(iterator seqIter,
				extDirection dir, uint8_t base);
		bool removeExtensionByIter(iterator seqIter,
				extDirection dir, SeqExt ext);

		/** Call the observers of the specified sequence. */
		void notify(const PackedSeq& seq)
		{
			if (m_seqObserver != NULL)
				m_seqObserver(this, seq);
		}

		/** The underlying collection. */
		SequenceDataHash* m_pSequences;

		/** The observers. Only a single observer is implemented.*/
		SeqObserver m_seqObserver;

		/** Whether adjacency information has been loaded. */
		bool m_adjacencyLoaded;
};

#endif
