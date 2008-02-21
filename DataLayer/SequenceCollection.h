#ifndef SequenceCollection_H
#define SequenceCollection_H

#include <deque>
#include <vector>
#include <set>
#include <stdio.h>
#include "ISequenceCollection.h"
#include "Sequence.h"
#include "PackedSeq.h"
#include "HitRecord.h"
#include "CommonDefs.h"

using namespace std;

enum SpaceState
{
	SS_LOADING,
	SS_FINALIZED,
	SS_READY
};


typedef std::deque<PackedSeq> SequenceData;
typedef SequenceData::iterator SequenceCollectionIter;
typedef SequenceData::iterator ConstSequenceCollectionIter;

typedef std::pair<SequenceCollectionIter, SequenceCollectionIter> SequenceIterPair;


// This class implements a collection of PackedSeqs with functions to manipulate the data
// It is meant to be a storage class only and should have minimal logic
class SequenceCollection : public ISequenceCollection
{
	public:
	
		//Allocates phase space
		SequenceCollection();
		
		//Deallocates phase space
		~SequenceCollection();

		// add a single sequence to the collection
		void addSequence(const PackedSeq& seq);
		
		// remove a sequence from the collection
		void removeSequence(const PackedSeq& seq);
				
		// end the data load and make the sequence space ready for data read
		void finalize();
		
		// check if a sequence exists
		bool checkForSequence(const PackedSeq& seq) const;
		
		// Set flag for sequence seq
		void markSequence(const PackedSeq& seq, SeqFlag flag);
		
		// Find if this sequence has the specified flag set
		bool checkSequenceFlag(const PackedSeq& seq, SeqFlag flag);
		
		// remove the extension to the sequence
		void removeExtension(const PackedSeq& seq, extDirection dir, char base);
	
		// Get the iterator pointing to the first sequence in the bin
		SequenceCollectionIter getStartIter() const;
		
		// Get the iterator pointing to the last sequence in the bin
		SequenceCollectionIter getEndIter() const;
		
		// does this sequence extend from a different node?
		bool hasParent(const PackedSeq& seq) const;

		// does this sequence have an extension?
		bool hasChild(const PackedSeq& seq) const;
		
		// Return the number of sequences in the collection
		int countAll() const;
						

	private:

		// Functions to get iterators to the sequence
		
		// Get the iterator to the sequence and its reverse complement
		// If they don't exist m_pSequences->end() will be returned in the iterator
		SequenceIterPair GetSequenceIterators(const PackedSeq& seq) const;
		SequenceCollectionIter FindSequence(const PackedSeq& seq) const;
		
		// Check if duplicate entries exist
		bool checkForDuplicates() const;		
		
		// Iterator versions of get/set functions
		// These should only be called from this class, hence they are private
		void removeSequencePrivate(SequenceIterPair seqIters);
		bool hasChildPrivate(SequenceCollectionIter seqIter) const;
		bool hasParentPrivate(SequenceCollectionIter seqIter) const;
		bool checkSequenceFlagPrivate(SequenceCollectionIter& seqIter, SeqFlag flag);
		void markSequencePrivate(SequenceCollectionIter& seqIter, SeqFlag flag);
		void removeExtensionPrivate(SequenceCollectionIter& seqIter, extDirection dir, char base);
		bool checkSequenceExtensionPrivate(SequenceCollectionIter& seqIter, extDirection dir, char base) const;
		
		// Data members
		
		// pointer to the actual collection (typedef'd above)
		SequenceData* m_pSequences;
		
		// the state of the space
		SpaceState m_state;
};

#endif
