#ifndef ISEQUENCECOLLECTION_H
#define ISEQUENCECOLLECTION_H

#include "CommonDefs.h"
#include "CommonUtils.h"
#include "PackedSeq.h"
#include "HitRecord.h"

// Interface class for a sequence space (collection of sequences)
// This pure virtual class defines the minimum set of functions a sequence space must provide

class ISequenceCollection
{
	public:
		virtual ~ISequenceCollection() {};
		
		// add a sequence to the collection
		virtual void addSequence(const PackedSeq& seq) = 0;
		
		// remove a sequence from the collection
		virtual void removeSequence(const PackedSeq& seq) = 0;
		
		// end the data load and make the sequence space ready for data read
		virtual void finalize() = 0;
		
		// check if a sequence exists
		virtual bool checkForSequence(const PackedSeq& seq) const = 0;
		
		// Set flag for sequence seq
		virtual void markSequence(const PackedSeq& seq, SeqFlag flag) = 0;
		
		// Find if this sequence has the specified flag set
		virtual bool checkSequenceFlag(const PackedSeq& seq, SeqFlag flag) = 0;
		
		// does this sequence extend from a different node?
		virtual bool hasParent(const PackedSeq& seq) const = 0;

		// does this sequence have an extension?
		virtual bool hasChild(const PackedSeq& seq) const = 0;
		
		// Return the number of sequences in the collection
		virtual int countAll() const = 0;								
};

#endif
