#ifndef PAIRRECORD_H
#define PAIRRECORD_H

#include <map>
#include <vector>
#include "Sequence.h"

class PairRecord
{
	
	public:
		PairRecord(const std::vector<Sequence>& allseqs);
		
		bool checkForPairs(const Sequence& seq) const;
		
		const std::vector<Sequence>& getPairs(const Sequence& seq) const;
	
	private:
	
		// the pair lookup table, a hash of sequences pointing to a vector of sequences, their pairs
		std::map<Sequence, std::vector<Sequence> > m_pairLookup;
};

#endif
