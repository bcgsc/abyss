#ifndef _ROLLIINGBLOOMDBGVERTEX_SET_H_
#define _ROLLIINGBLOOMDBGVERTEX_SET_H_

#include "BloomDBG/RollingBloomDBG.h"
#include "Graph/Path.h"
#include "Common/Sequence.h"

#include<unordered_set>

class RollingBloomDBGVertexSet
{
protected:

	std::unordered_set<RollingBloomDBGVertex> m_set;
	unsigned m_k;

public:

	RollingBloomDBGVertexSet(unsigned k) : m_k(k) { }

	bool containsRollingBloomDBGVertex(const RollingBloomDBGVertex& RollingBloomDBGVertex)
	{
		return m_set.find(RollingBloomDBGVertex) != m_set.end();
	}

	void addRollingBloomDBGVertex(const RollingBloomDBGVertex& RollingBloomDBGVertex)
	{
		m_set.insert(RollingBloomDBGVertex);
	}

	void loadSeq(Sequence seq)
	{
		if (seq.length() < m_k)
			return;

		flattenAmbiguityCodes(seq);
		for (RollingHashIterator it(seq, 1, m_k); it != RollingHashIterator::end(); ++it) {
			RollingBloomDBGVertex curr(it.kmer().c_str(), it.rollingHash());
			m_set.insert(curr);
		}
	}
};

#endif