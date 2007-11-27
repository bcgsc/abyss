#ifndef HITRECORD_H
#define HITRECORD_H

#include <vector>
#include "Sequence.h"

class HitRecord
{
	public:
		HitRecord();
		
		// add a hit to the record
		void addHit(const Sequence& seq);
		
		// get the number of hits
		int getNumHits() const;
		
		// get a specific hit
		Sequence getHit(int num) const;
		
		// get the first  hit
		Sequence getFirstHit() const;
		
	private:
		std::vector<Sequence> m_hits;
	
};

#endif
