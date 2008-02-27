#ifndef ISEQUENCECOLLECTION_H
#define ISEQUENCECOLLECTION_H

#include "CommonDefs.h"
#include "CommonUtils.h"
#include "PackedSeq.h"
#include "HitRecord.h"

// Interface class for a sequence collection (the lowest level of storage of a large number of sequences)
// This pure virtual class defines the minimum set of functions a sequence collection must provide


// Most operations are performed on the forward and reverse reads simulatenously, this structure holds the result of such operations
struct ResultPair
{
	bool forward;
	bool reverse;
};

class ISequenceCollection
{
	public:
		virtual ~ISequenceCollection() {};
		
		// add a sequence to the collection
		virtual void add(const PackedSeq& seq) = 0;
		
		// remove a sequence from the collection
		virtual void remove(const PackedSeq& seq) = 0;
		
		// end the data load and make the sequence space ready for data read
		virtual void finalize() = 0;
		
		// check if a sequence exists
		virtual ResultPair exists(const PackedSeq& seq) const = 0;
		
		// Set flag for sequence seq
		virtual void setFlag(const PackedSeq& seq, SeqFlag flag) = 0;
		
		// Find if this sequence has the specified flag set
		virtual ResultPair checkFlag(const PackedSeq& seq, SeqFlag flag) = 0;
		
		// does this sequence extend from a different node?
		virtual bool hasParent(const PackedSeq& seq) const = 0;

		// does this sequence have an extension?
		virtual bool hasChild(const PackedSeq& seq) const = 0;
		
		// Return the number of sequences in the collection
		virtual int count() const = 0;
		
		// remove the extension to the sequence
		virtual void removeExtension(const PackedSeq& seq, extDirection dir, char base) = 0;
		
		// add an extension to the sequence
		virtual void setExtension(const PackedSeq& seq, extDirection dir, SeqExt extension) = 0;
		
		// check if the extension exists
		virtual ResultPair checkExtension(const PackedSeq& seq, extDirection dir, char base) const = 0;
};

#endif
