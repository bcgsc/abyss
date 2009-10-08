#ifndef SEQUENCECOLLECTIONHASH_H
#define SEQUENCECOLLECTIONHASH_H 1

#include "ISequenceCollection.h"
#include "PackedSeq.h"

typedef std::pair<bool, SeqExt> SeqExtResult;

// This class implements a collection of PackedSeqs with functions to manipulate the data
// It is meant to be a storage class only and should have minimal logic for manipulating the data except for getters/setters
class SequenceCollectionHash : public ISequenceCollection
{
	public:
	
		//Allocates phase space
		SequenceCollectionHash();
		
		//Deallocates phase space
		~SequenceCollectionHash();

		// add a single sequence to the collection
		void add(const PackedSeq& seq);

		/** Remove the specified sequence if it exists. */
		void remove(const PackedSeq& seq)
		{
			setFlag(seq, SF_DELETE);
		}

		// Clean up by erasing sequences flagged as deleted.
		unsigned cleanup();

		// Print the load of the hash table.
		void printLoad() const;

		// Set flag for sequence seq
		void setFlag(const PackedSeq& seq, SeqFlag flag);

		// Clear the specified flag from every sequence in the collection
		void wipeFlag(SeqFlag flag);

		// set a base extension
		bool setBaseExtension(const PackedSeq& seq, extDirection dir,
				uint8_t base);

		// remove the extension to the sequence
		bool removeExtension(const PackedSeq& seq,
				extDirection dir, uint8_t base);

		// clear the extensions of the sequence
		void clearExtensions(const PackedSeq& seq, extDirection dir);

		// get the extensions of a sequence
		bool getSeqData(const PackedSeq& seq,
				ExtensionRecord& extRecord, int& multiplicity) const;

		const PackedSeq& getSeqAndData(const PackedSeq& key) const;

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
		SequenceHashIterPair GetSequenceIterators(const PackedSeq& seq) const;
		const_iterator SequenceCollectionHash::find(
				const PackedSeq& key, bool& rc) const;
		const PackedSeq& getSeqAndData(
				const SequenceHashIterPair& iters) const;

		bool setBaseExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, uint8_t base);
		bool removeExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, uint8_t base);
		bool clearExtensionsByIter(SequenceCollectionHashIter& seqIter, extDirection dir);

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
