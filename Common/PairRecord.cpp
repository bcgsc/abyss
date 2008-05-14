#include "PairRecord.h"

PairRecord::PairRecord()
{

}

void PairRecord::addPairs(const PackedSeq& seq1, const PackedSeq& seq2)
{
	m_pairLookup[seq1].push_back(seq2);
	m_pairLookup[seq2].push_back(seq1);	
}

bool PairRecord::checkForPairs(const PackedSeq& seq) const
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

const PSequenceVector PairRecord::getPairs(const PackedSeq& seq) const
{
	PSequenceVector ret;
	addPairs(seq, ret);
	addPairs(reverseComplement(seq), ret);
	return ret;
}

void PairRecord::addPairs(const PackedSeq& seq, PSequenceVector& vec) const
{
	std::map<PackedSeq, PSequenceVector >::const_iterator findIter = m_pairLookup.find(seq);
	if(findIter != m_pairLookup.end())
	{
		for(PSequenceVector::const_iterator iter = findIter->second.begin(); iter != findIter->second.end(); ++iter)
		{
			vec.push_back(*iter);
		}
	}
	return;
}
