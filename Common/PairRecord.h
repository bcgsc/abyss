#ifndef PAIRRECORD_H
#define PAIRRECORD_H

#include <map>
#include <vector>
#include "PackedSeq.h"

class PairRecord
{
	
	public:
		PairRecord(const PSequenceVector& allseqs);
		
		bool checkForPairs(const PackedSeq& seq) const;
		
		const PSequenceVector& getPairs(const PackedSeq& seq) const;
	
	private:
	
		// the pair lookup table, a hash of sequences pointing to a vector of sequences, their pairs
		std::map<PackedSeq, PSequenceVector > m_pairLookup;
};

#endif
