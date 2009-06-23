#ifndef HITRECORD_H
#define HITRECORD_H

#include "PackedSeq.h"
#include <vector>

typedef std::vector<PackedSeq> HitVector;

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
		HitVector m_hits;
	
};

#endif
