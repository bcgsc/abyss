#ifndef SEQRECORD_H
#define SEQRECORD_H

#include <map>
#include "CommonDefs.h"
#include "Sequence.h"

// Class to record which sequences have been seen and extended (to prevent duplicate extensions)
class SeqRecord
{
	public:
		SeqRecord();
		
		// add a lot of sequences to the record via a vector
		void addMultipleSequences(const std::vector<Sequence>& seqVec);
		
		// add a lot of sequences to the record
		void addMultipleSequences(const SequenceMap& seqMap);
		
		// add a single sequence to the record
		void addSequence(const Sequence& seq);
		
		// get the number of times this sequence appears in the record
		int getMultiplicity(const Sequence& seq) const;
		
		// check if this sequence is contained in the record
		bool contains(const Sequence& seq) const;
		
	private:
		std::map<Sequence, int> m_seqRecord;
};

#endif
