#ifndef SEQUENCECOLLECTION_H
#define SEQUENCECOLLECTION_H

#include <deque>
#include <vector>
#include <set>
#include <stdio.h>
#include "ISequenceCollection.h"
#include "Sequence.h"
#include "PackedSeq.h"
#include "HitRecord.h"
#include "CommonDefs.h"

typedef std::vector<PackedSeq> SequenceData;
typedef SequenceData::iterator SequenceCollectionIter;
typedef SequenceData::iterator ConstSequenceCollectionIter;

typedef std::pair<SequenceCollectionIter, SequenceCollectionIter> SequenceIterPair;


// This class implements a collection of PackedSeqs with functions to manipulate the data
// It is meant to be a storage class only and should have minimal logic for manipulating the data except for getters/setters
class SequenceCollection : public ISequenceCollection
{
	public:
	
		//Allocates phase space
		SequenceCollection();
		
		//Deallocates phase space
		~SequenceCollection();

		// add a single sequence to the collection
		void add(const PackedSeq& seq);
		
		// remove a sequence from the collection
		void remove(const PackedSeq& seq);
				
		// end the data load and make the sequence space ready for data read
		void finalize();
		
		// check if a sequence exists
		bool exists(const PackedSeq& seq);
		
		// Set flag for sequence seq
		void setFlag(const PackedSeq& seq, SeqFlag flag);
		
		// Find if this sequence has the specified flag set
		bool checkFlag(const PackedSeq& seq, SeqFlag flag);
		
		// add extension
		void setExtension(const PackedSeq& seq, extDirection dir, SeqExt extension);
		
		// remove the extension to the sequence
		void removeExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// check if the extension exists
		bool checkExtension(const PackedSeq& seq, extDirection dir, char base);
	
		// Get the iterator pointing to the first sequence in the bin
		SequenceCollectionIter getStartIter();
		
		// Get the iterator pointing to the last sequence in the bin
		SequenceCollectionIter getEndIter();
		
		// does this sequence extend from a different node?
		bool hasParent(const PackedSeq& seq);

		// does this sequence have an extension?
		bool hasChild(const PackedSeq& seq);
		
		// Return the number of sequences in the collection
		int count() const;
		
		// Dummy function to conform to the interface
		APResult pumpNetwork() { return APR_NONE; }
						

	private:

		// Functions to get iterators to the sequence
		
		// Get the iterator to the sequence and its reverse complement
		// If they don't exist m_pSequences->end() will be returned in the iterator
		SequenceIterPair GetSequenceIterators(const PackedSeq& seq) const;
		SequenceCollectionIter FindSequence(const PackedSeq& seq) const;
		
		// Check if duplicate entries exist
		bool checkForDuplicates() const;		
		
		// Iterator versions of modification functions
		// These should only be called from this class, hence they are private
		void removeByIter(SequenceIterPair seqIters);
		bool hasChildByIter(SequenceCollectionIter seqIter) const;
		bool hasParentByIter(SequenceCollectionIter seqIter) const;
		bool checkFlagByIter(SequenceCollectionIter& seqIter, SeqFlag flag);
		void setFlagByIter(SequenceCollectionIter& seqIter, SeqFlag flag);
		void setExtensionByIter(SequenceCollectionIter& seqIter, extDirection dir, SeqExt extension);
		void removeExtensionByIter(SequenceCollectionIter& seqIter, extDirection dir, char base);
		bool checkExtensionByIter(SequenceCollectionIter& seqIter, extDirection dir, char base) const;
		bool existsByIter(SequenceCollectionIter& seqIter) const;
		
		
		// Data members
		
		// pointer to the actual collection (typedef'd above)
		SequenceData* m_pSequences;
		
		// the state of the space
		CollectionState m_state;
};

#endif
