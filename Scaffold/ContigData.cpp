#include "ContigData.h"

ContigData::ContigData(const Sequence& s, size_t kmer) : seq(s), m_kmer(kmer)
{
	// Create the seq set
	for(size_t idx = 0; idx < seq.length() - m_kmer; ++idx)
	{
		PackedSeq insSeq(seq.substr(idx, m_kmer));
		m_seqSet.insert(insSeq);	
		m_seqSet.insert(reverseComplement(insSeq));
	}
}
