#ifndef HITRECORD_H
#define HITRECORD_H

#include "CommonDefs.h"
#include "PackedSeq.h"

class HitRecord
{
	public:
		HitRecord();
		
		// add a hit to the record
		void addHit(const PackedSeq& seq);
		
		// get the number of hits
		int getNumHits() const;
		
		// get a specific hit
		PackedSeq getHit(int num) const;
		
		// get the first  hit
		PackedSeq getFirstHit() const;
		
	private:
		PSequenceVector m_hits;
	
};

#endif
