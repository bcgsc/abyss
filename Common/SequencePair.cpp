#include "SequencePair.h"

SequencePair::SequencePair(Sequence s, int pos) : m_seq(s), m_maxPos(pos), m_distance(0), m_stddev(0.0f)
{
	
}

Sequence SequencePair::getSequence() const
{
	return m_seq;
}

int SequencePair::getMaxPos() const
{
	return m_maxPos;	
}
	
int SequencePair::getDistance() const
{
	return m_distance;	
}
		
double SequencePair::getStdDev() const
{
	return m_stddev;	
}
