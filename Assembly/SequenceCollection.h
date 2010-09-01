#ifndef SEQUENCECOLLECTION_H
#define SEQUENCECOLLECTION_H 1

#include "ISequenceCollection.h"
#include <cassert>

/** A map of Kmer to KmerData. */
class SequenceCollectionHash : public ISequenceCollection
{
	public:
		SequenceCollectionHash();

		void add(const Kmer& seq);

		/** Remove the specified sequence if it exists. */
		void remove(const Kmer& seq)
		{
			setFlag(seq, SF_DELETE);
		}

		// Clean up by erasing sequences flagged as deleted.
		unsigned cleanup();

		/** Shrink the hash table. */
		void shrink() {
#if USING_EXT_HASH_MAP
			m_data.resize(0);
#else
			m_data.rehash(0);
#endif
			printLoad();
		}

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

		const value_type& getSeqAndData(const Kmer& key) const;

		iterator begin() { return m_data.begin(); }
		const_iterator begin() const { return m_data.begin(); }
		iterator end() { return m_data.end(); }
		const_iterator end() const { return m_data.end(); }

		/** Return true if this collection is empty. */
		bool empty() const { return m_data.empty(); }

		/** Return the number of sequences in this collection. */
		size_t size() const { return m_data.size(); }

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
			(void)f;
			m_seqObserver = NULL;
		}

		void load(const char *path);
		void store(const char* path);
		bool isAdjacencyLoaded() const { return m_adjacencyLoaded; }
		void setColourSpace(bool flag);

	private:
		iterator find(const Kmer& key) { return m_data.find(key); }
		const_iterator find(const Kmer& key) const
		{
			return m_data.find(key);
		}

		iterator find(const Kmer& key, bool& rc);
		const_iterator find(const Kmer& key, bool& rc) const;

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

#endif
