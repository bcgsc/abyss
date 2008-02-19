#ifndef SIMPLESEQUENCESPACE_H
#define SIMPLESEQUENCESPACE_H

#include <deque>
#include <vector>
#include <set>
#include <stdio.h>
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


typedef std::deque<PackedSeq> SequenceCollection;
typedef SequenceCollection::iterator SequenceCollectionIter;
typedef SequenceCollection::iterator ConstSequenceCollectionIter;

typedef std::pair<SequenceCollectionIter, SequenceCollectionIter> SequenceIterPair;

class SimpleSequenceSpace
{
	public:
	
		//Allocates phase space
		SimpleSequenceSpace(int readLength, Coord4 startCoord, Coord4 size);
		
		//Deallocates phase space
		~SimpleSequenceSpace();

		// add a single sequence
		void addSequence(const PackedSeq& seq, bool boundsCheck = false);
		
		// remove a sequence, the PackedSeq version is O(logn), the iterator version is O(n)
		void removeSequence(const PackedSeq& seq);
		
		// remove the extension to the sequence
		void removeExtension(const PackedSeq& seq, extDirection dir, char base);
		
		// trim and sort the vectors
		void finalizeBins();
		
		// generate the adjacency info for every sequence in the collection
		void generateAdjacency();
		
		// Check if duplicate entries exist
		bool checkForDuplicates() const;
		
		// get the multiplicity of the sequence
		int getMultiplicity(const PackedSeq& seq);
		
		// check if a sequence exists
		bool checkForSequence(const PackedSeq& seq) const;
		
		// Mark a sequence in the phase space, the PackedSeq version is O(logn), the iterator version is O(n)
		void markSequence(const PackedSeq& seq, SeqFlag flag);
		
		// Find if this sequence has the specified flag set, the PackedSeq version is O(logn), the iterator version is O(n)
		bool checkSequenceFlag(const PackedSeq& seq, SeqFlag flag);

		// calculate whether this sequence has an extension in the phase space
		HitRecord calculateExtension(const PackedSeq& currSeq, extDirection dir) const;
		
		// Get the iterator pointing to the first sequence in the bin
		SequenceCollectionIter getStartIter() const;
		
		// Get the iterator pointing to the last sequence in the bin
		SequenceCollectionIter getEndIter() const;
		
		// does this sequence extend from a different node?
		bool hasParent(const PackedSeq& seq) const;

		// does this sequence have an extension?
		bool hasChild(const PackedSeq& seq) const;

		
		// print everything
		void printAll() const;
		
		int countAll() const;
						

	private:

		SequenceIterPair GetSequenceIterators(const PackedSeq& seq) const;
		SequenceCollectionIter FindSequence(const PackedSeq& seq) const;
		
		// Iterator versions of get/set functions
		// These are not exposed to the public
		void removeSequencePrivate(SequenceIterPair seqIters);
		bool hasChildPrivate(SequenceCollectionIter seqIter) const;
		bool hasParentPrivate(SequenceCollectionIter seqIter) const;
		bool checkSequenceFlagPrivate(SequenceCollectionIter& seqIter, SeqFlag flag);
		void markSequencePrivate(SequenceCollectionIter& seqIter, SeqFlag flag);
		void removeExtensionPrivate(SequenceCollectionIter& seqIter, extDirection dir, char base);
		bool checkSequenceExtensionPrivate(SequenceCollectionIter& seqIter, extDirection dir, char base) const;
		
		SimpleSequenceSpace();
		
		SequenceCollection* m_pSequences;
		Coord4 m_minCoord;
		Coord4 m_maxCoord;
		Coord4 m_size;
		int m_readLength;
		SpaceState m_state;
};

#endif
