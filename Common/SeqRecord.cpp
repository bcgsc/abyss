#include "SeqRecord.h"

// Record of the multiplicity (and hence existence of sequences)

SeqRecord::SeqRecord()
{
	
}

// add a collection of sequences in a vector
void SeqRecord::addMultipleSequences(const PSequenceVector& seqVec)
{
	for(ConstPSequenceVectorIterator iter = seqVec.begin(); iter != seqVec.end(); iter++)
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
void SeqRecord::addSequence(const PackedSeq& seq)
{
	m_seqRecord[seq]++;
	PackedSeq reverseComp(seq);
	reverseComp.reverseComplement();
	m_seqRecord[reverseComp]++;	
}


// return the multiplicity of a sequence	
int SeqRecord::getMultiplicity(const PackedSeq& seq) const
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
bool SeqRecord::contains(const PackedSeq& seq) const
{
	return getMultiplicity(seq) > 0;
}
