#include "SeqRecord.h"

// Record of the multiplicity (and hence existence of sequences)

SeqRecord::SeqRecord()
{
	
}

// add a collection of sequences in a vector
void SeqRecord::addMultipleSequences(const std::vector<Sequence>& seqVec)
{
	for(const_seq_iter iter = seqVec.begin(); iter != seqVec.end(); iter++)
	{
		addSequence(*iter);
	}	
}

// add a collection of sequences in a map
void SeqRecord::addMultipleSequences(const SequenceMap& seqMap)
{
	for(SequenceMap::const_iterator iter = seqMap.begin(); iter != seqMap.end(); iter++)
	{
		addSequence(iter->second);
	}	
}

// add a single sequence
void SeqRecord::addSequence(const Sequence& seq)
{
	m_seqRecord[seq]++;
	m_seqRecord[reverseComplement(seq)]++;	
}


// return the multiplicity of a sequence	
int SeqRecord::getMultiplicity(const Sequence& seq) const
{
	ConstSeqRecordIter iter = m_seqRecord.find(seq);
	if(iter != m_seqRecord.end())
	{
		return iter->second;
	}
	else
	{
		return 0;
	}	
}

// check if this sequence exists
bool SeqRecord::contains(const Sequence& seq) const
{
	return getMultiplicity(seq) > 0;
}
