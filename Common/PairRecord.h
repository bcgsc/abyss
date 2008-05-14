#ifndef PAIRRECORD_H
#define PAIRRECORD_H

#include <map>
#include <vector>
#include "PackedSeq.h"

class PairRecord
{
	
	public:
		PairRecord();
		
		void addPairs(const PackedSeq& seq1, const PackedSeq& seq2);
		bool checkForPairs(const PackedSeq& seq) const;
		
		const PSequenceVector getPairs(const PackedSeq& seq) const;
	
	private:
	
		void addPairs(const PackedSeq& seq, PSequenceVector& vec) const;
		
		// the pair lookup table, a hash of sequences pointing to a vector of sequences, their pairs
		std::map<PackedSeq, PSequenceVector > m_pairLookup;
};

#endif
