#ifndef SEQUENCECOLLECTIONHASH_H
#define SEQUENCECOLLECTIONHASH_H

#include <deque>
#include <vector>
#include <set>
#include <stdio.h>
#include <ext/hash_set>
#include "ISequenceCollection.h"
#include "Sequence.h"
#include "PackedSeq.h"
#include "HitRecord.h"
#include "CommonDefs.h"

using namespace std;

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
				
		// Get the multiplicity of a sequence
		int getMultiplicity(const PackedSeq& seq);
		
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
		
		// set a base extension
		void setBaseExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// remove the extension to the sequence
		void removeExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// clear the extensions of the sequence
		void clearExtensions(const PackedSeq& seq, extDirection dir);
		
		// check if the extension exists
		ResultPair checkExtension(const PackedSeq& seq, extDirection dir, char base);

		// Get the iterator pointing to the first sequence in the bin
		SequenceCollectionHashIter getStartIter() const;
		
		// Get the iterator pointing to the last sequence in the bin
		SequenceCollectionHashIter getEndIter() const;
		
		// does this sequence extend from a different node?
		bool hasParent(const PackedSeq& seq);

		// does this sequence have an extension?
		bool hasChild(const PackedSeq& seq);
		
		// Generate the initial cache of branch ends
		void cacheBranchEnds();
		
		// Make a copy of the branch end cache for the higher-level algorithms to operate on
		void copyBranchCache(PSeqSet& outset);
		
		// Return the number of sequences in the collection
		int count() const;
		
		// Pump the network. For this sequence collection (non-network) this function just returns
		virtual APResult pumpNetwork() { return APR_NONE; }
								
	private:


		// Functions to get iterators to the sequence
		
		// Get the iterator to the sequence and its reverse complement
		// If they don't exist m_pSequences->end() will be returned in the iterator
		SequenceHashIterPair GetSequenceIterators(const PackedSeq& seq) const;
		SequenceCollectionHashIter FindSequence(const PackedSeq& seq) const;
		
		// Check if duplicate entries exist
		bool checkForDuplicates() const;		
		
		// Iterator versions of modification functions
		// These should only be called from this class, hence they are private
		void removeByIter(SequenceHashIterPair seqIters);
		
		bool hasChildByIter(SequenceCollectionHashIter seqIter) const;
		bool hasParentByIter(SequenceCollectionHashIter seqIter) const;
		bool checkFlagByIter(SequenceCollectionHashIter& seqIter, SeqFlag flag);
		void setFlagByIter(SequenceCollectionHashIter& seqIter, SeqFlag flag);
		void setExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, SeqExt extension);
		void setBaseExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base);
		void removeExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base);
		void clearExtensionsByIter(SequenceCollectionHashIter& seqIter, extDirection dir);
		bool checkExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base) const;
		bool existsByIter(SequenceCollectionHashIter& seqIter) const;

		// Data members
		
		// pointer to the actual collection (typedef'd above)
		SequenceDataHash* m_pSequences;
		
		// the state of the space
		CollectionState m_state;
		
		// Cache of branch ends (sequences that have no extension on either side)
		PSeqSet m_branchEndCache;
};

#endif
