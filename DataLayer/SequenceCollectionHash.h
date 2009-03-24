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
		
		// remove a sequence from the collection
		void remove(const PackedSeq& seq);

		// Remove marked sequences from the collection.
		unsigned removeMarked();

		// Get the multiplicity of a sequence
		int getMultiplicity(const PackedSeq& seq);
		
		// Print the load of the hash table.
		void printLoad();
		
		// check if a sequence exists
		bool exists(const PackedSeq& seq);
		
		// Set flag for sequence seq
		void setFlag(const PackedSeq& seq, SeqFlag flag);
		
		// Find if this sequence has the specified flag set
		bool checkFlag(const PackedSeq& seq, SeqFlag flag);
		
		// Clear the specified flag from every sequence in the collection
		void wipeFlag(SeqFlag flag);

		// set a base extension
		bool setBaseExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// remove the extension to the sequence
		bool removeExtension(const PackedSeq& seq,
				extDirection dir, char base);
		
		// clear the extensions of the sequence
		void clearExtensions(const PackedSeq& seq, extDirection dir);
		
		// check if the extension exists
		ResultPair checkExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// get the extensions of a sequence
		bool getSeqData(const PackedSeq& seq, ExtensionRecord& extRecord, int& multiplicity);

		const PackedSeq& getSeqAndData(const PackedSeq& key) const;

		// Get the iterator pointing to the first sequence in the bin
		SequenceCollectionHashIter getStartIter() const;
		
		// Get the iterator pointing to the last sequence in the bin
		SequenceCollectionHashIter getEndIter() const;
		
		// does this sequence extend from a different node?
		bool hasParent(const PackedSeq& seq);

		// does this sequence have an extension?
		bool hasChild(const PackedSeq& seq);
		
		// Return the number of sequences in the collection
		size_t count() const;

		// Not a network sequence collection. Nothing to do.
		virtual unsigned pumpNetwork() { return 0; }

		/** Attach the specified observer. */
		virtual void attach(SeqObserver f)
		{
			assert(m_seqObserver == NULL);
			m_seqObserver = f;
		}

		/** Detach the specified observer. */
		virtual void detach(SeqObserver f)
		{
			assert(m_seqObserver == f);
			m_seqObserver = NULL;
		}

	private:
		// Functions to get iterators to the sequence
		
		// Get the iterator to the sequence and its reverse complement
		// If they don't exist m_pSequences->end() will be returned in the iterator
		SequenceHashIterPair GetSequenceIterators(const PackedSeq& seq) const;
		SequenceCollectionHashIter FindSequence(const PackedSeq& seq) const;
		const PackedSeq& getSeqAndData(
				const SequenceHashIterPair& iters) const;

		// Iterator versions of modification functions
		// These should only be called from this class, hence they are private
		void removeByIter(SequenceHashIterPair seqIters);
		
		bool hasChildByIter(SequenceCollectionHashIter seqIter) const;
		bool hasParentByIter(SequenceCollectionHashIter seqIter) const;
		bool checkFlagByIter(SequenceCollectionHashIter& seqIter, SeqFlag flag);
		void setFlagByIter(SequenceCollectionHashIter& seqIter, SeqFlag flag);
		bool setBaseExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base);
		void removeExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base);
		void clearExtensionsByIter(SequenceCollectionHashIter& seqIter, extDirection dir);
		bool checkExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base) const;
		bool existsByIter(SequenceCollectionHashIter& seqIter) const;
		SeqExt getExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir) const;

		/** Call the observers of the specified sequence. */
		void notify(const PackedSeq& seq)
		{
			if (m_seqObserver != NULL)
				m_seqObserver(this, seq);
		}

		// Data members
		
		// pointer to the actual collection (typedef'd above)
		SequenceDataHash* m_pSequences;

		/** The observers. Only a single observer is implemented.*/
		SeqObserver m_seqObserver;
};

#endif
