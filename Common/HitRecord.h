#ifndef HITRECORD_H
#define HITRECORD_H

#include "CommonDefs.h"
#include "PackedSeq.h"


struct Hit
{
	Hit(PackedSeq s, bool b) : seq(s), isReverse(b) {}
	PackedSeq seq;
	bool isReverse;
};

typedef std::vector<Hit> HitVector;

class HitRecord
{
	public:
		HitRecord();
		
		// add a hit to the record
		void addHit(const PackedSeq& seq, bool isReverse);
		
		// get the number of hits
		int getNumHits() const;
		
		// get a specific hit
		Hit getHit(int num) const;
		
		// get the first  hit
		Hit getFirstHit() const;
		
	private:
		HitVector m_hits;
	
};

#endif
