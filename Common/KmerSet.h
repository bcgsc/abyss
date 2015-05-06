#ifndef _KMER_SET_H_
#define _KMER_SET_H_

#include "Common/Kmer.h"
#include "Graph/Path.h"
#include "Common/Sequence.h"

class KmerSet
{
protected:

	unordered_set<Kmer> m_set;
	unsigned m_k;

public:

	KmerSet(unsigned k) : m_k(k) { }

	bool containsKmer(const Kmer& kmer)
	{
		return m_set.find(kmer) != m_set.end();
	}

	void addKmer(const Kmer& kmer)
	{
		m_set.insert(kmer);
	}

	void loadSeq(Sequence seq)
	{
		if (seq.length() < m_k)
			return;

		flattenAmbiguityCodes(seq);
		for (unsigned i = 0; i < seq.length() - m_k + 1; ++i) {
			std::string kmerStr = seq.substr(i, m_k);
			size_t pos = kmerStr.find_first_not_of("AGCTagct");
			if (pos != std::string::npos) {
				i += pos;
				continue;
			}
			m_set.insert(Kmer(kmerStr));
		}
	}
};

#endif
