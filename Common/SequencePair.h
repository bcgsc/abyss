#ifndef SEQUENCEPAIR_H
#define SEQUENCEPAIR_H

#include "PackedSeq.h"

class SequencePair
{
	public:
	
		SequencePair(PackedSeq s, int pos);
			
		// get the sequence
		PackedSeq getSequence() const;
		
		// get the maximum position that we can expect this pair
		int getMaxPos() const;
		
		// get the expected distance from its pair
		int getDistance() const;
		
		// get the standard deviation
		double getStdDev() const;
		
	private:
		PackedSeq m_seq;
		int m_maxPos;
		int m_distance;
		double m_stddev;
};

#endif
