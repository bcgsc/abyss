#ifndef SEQRECORD_H
#define SEQRECORD_H

#include <map>
#include "PackedSeq.h"
#include "Sequence.h"

typedef std::map<PackedSeq, char>::const_iterator ConstSeqRecordIter;

// Class to record which sequences have been seen and extended (to prevent duplicate extensions)
class SeqRecord
{
	public:
		SeqRecord();
		
		// add a lot of sequences to the record via a vector
		void addMultipleSequences(const PSequenceVector& seqVec);
		
		// add a lot of sequences to the record
		void addMultipleSequences(const SequenceMap& seqMap);
		
		// add a single sequence to the record
		void addSequence(const PackedSeq& seq);
		
		// get the number of times this sequence appears in the record
		int getMultiplicity(const PackedSeq& seq) const;
		
		// check if this sequence is contained in the record
		bool contains(const PackedSeq& seq) const;
		
		
		// Methods that pass through to the underlying stl type
		int size() const { return m_seqRecord.size(); }
		
		ConstSeqRecordIter begin() const { return m_seqRecord.begin(); }
		ConstSeqRecordIter end() const { return m_seqRecord.end(); }
				
	private:
	
		int getMultiplicityInternal(const PackedSeq& seq) const;
		
		std::map<PackedSeq, char> m_seqRecord;
};

#endif
