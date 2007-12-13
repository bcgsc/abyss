#ifndef SEQRECORD_H
#define SEQRECORD_H

#include <ext/hash_map>
#include "CommonDefs.h"
#include "Sequence.h"

using namespace std;
using namespace __gnu_cxx;

typedef hash_map<Sequence, int>::const_iterator ConstSeqRecordIter;

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
		
		
		// Methods that pass through to the underlying stl type
		int size() const { return m_seqRecord.size(); }
		
		ConstSeqRecordIter begin() const { return m_seqRecord.begin(); }
		ConstSeqRecordIter end() const { return m_seqRecord.end(); }
				
	private:
	
		int getMultiplicityInternal(const Sequence& seq) const;
		
		hash_map<Sequence, int> m_seqRecord;
};

#endif
