#include "HitRecord.h"

HitRecord::HitRecord()
{
	
	
}

// add a hit to the record
void HitRecord::addHit(const PackedSeq& seq, bool isReverse)
{
	m_hits.push_back(Hit(seq, isReverse));
}

// get the number of hits
int HitRecord::getNumHits() const
{
	return m_hits.size();
}

Hit HitRecord::getHit(int num) const
{
	return m_hits[num];
}

// get the first  hit
Hit HitRecord::getFirstHit() const
{
	return m_hits.front();
}
