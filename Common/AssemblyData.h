#ifndef ASSEMBLYDATA_H
#define ASSEMBLYDATA_H

#include <deque>
#include <vector>
#include <set>
#include <stdio.h>
#include "Sequence.h"
#include "PackedSeq.h"
#include "HitRecord.h"
#include "CommonDefs.h"
#include "SequenceCollection.h"

enum DataState
{
	DS_LOADING,
	DS_READY
};

class AssemblyData
{
	public:
	
		//Allocates phase space
		AssemblyData();
		
		//Deallocates phase space
		~AssemblyData();

		// add a single sequence to the collection
		void addSequence(const PackedSeq& seq);
		
		// remove a sequence from the collection
		void removeSequence(const PackedSeq& seq);
				
		// end the data load and make the sequence space ready for data read
		void finalize();
		
		// generate the adjacency info for every sequence in the collection
		void generateAdjacency();
		
		// check if a sequence exists
		bool checkForSequence(const PackedSeq& seq) const;
		
		// Set flag for sequence seq
		void markSequence(const PackedSeq& seq, SeqFlag flag);
		
		// Find if this sequence has the specified flag set
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
		
		// Return the number of sequences in the collection
		int countAll() const;
						

	private:

		// Data members
		
		// pointer to the actual collection (typedef'd above)
		SequenceCollection* m_pSequences;
		
		// The state of the data set
		// DS_LOADING = Don't read from the data set, load only
		// READY = ok to start reading and processing
		DataState m_state;
};

#endif
