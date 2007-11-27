#include "PairRecord.h"

PairRecord::PairRecord(const std::vector<Sequence>& allseqs)
{
	const_seq_iter pairIter;
	
	// the pairs are located sequentially in the vector
	for(const_seq_iter iter = allseqs.begin(); iter != allseqs.end(); iter += 2)
	{
		pairIter = iter + 1;
		m_pairLookup[*iter].push_back(*pairIter);
		m_pairLookup[*pairIter].push_back(*iter);
	}
}

bool PairRecord::checkForPairs(const Sequence& seq) const
{
	if(m_pairLookup.count(seq) > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

const std::vector<Sequence>& PairRecord::getPairs(const Sequence& seq) const
{
	// for any sequence, there HAS to be some pairs
	assert(m_pairLookup.count(seq) > 0);
	return m_pairLookup.find(seq)->second;
}
