#include "HitRecord.h"

HitRecord::HitRecord()
{
	
	
}

// add a hit to the record
void HitRecord::addHit(const PackedSeq& seq)
{
	m_hits.push_back(seq);
}

// get the number of hits
int HitRecord::getNumHits() const
{
	return m_hits.size();
}

PackedSeq HitRecord::getHit(int num) const
{
	return m_hits[num];
}

// get the first  hit
PackedSeq HitRecord::getFirstHit() const
{
	return m_hits.front();
}
