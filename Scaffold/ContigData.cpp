#include "ContigData.h"

ContigData::ContigData(const ContigID& id, const Sequence& s, size_t kmer, int copyNumber) : m_seq(s), m_kmer(kmer), m_copyNumber(copyNumber)
{
	// Create the seq set
	m_kmerVec.reserve(m_seq.length() - m_kmer + 1);
	
	// Add the initial sequences
	addSeqs(m_seq, m_kmerVec.begin(), true);
	
	addID(id, SENSE);
	addID(id, ANTISENSE);
}

void ContigData::addSeqs(const Sequence& s, CKDataVec::iterator position, bool usable)
{
	CKDataVec tempVec;
	for(size_t idx = 0; idx < s.length() - m_kmer + 1; ++idx)
	{
		PackedSeq insSeq(s.substr(idx, m_kmer));
		ContigKmerData ndata;
		ndata.seq = insSeq;
		ndata.usable = usable;
		tempVec.push_back(ndata);
	}	
	
	m_kmerVec.insert(position, tempVec.begin(), tempVec.end());
	
}

void ContigData::validateSequences() const
{
	assert((m_seq.length() - m_kmer + 1) == m_kmerVec.size());
	
	for(size_t idx = 0; idx < m_seq.length() - m_kmer + 1; ++idx)
	{
		PackedSeq refSeq(m_seq.substr(idx, m_kmer));
		const PackedSeq& kmer = m_kmerVec[idx].seq;
		
		assert(refSeq == kmer);
	}
}

void ContigData::copySeqSet(PSeqSet& outSeqs) const
{
	for(CKDataVec::const_iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		outSeqs.insert(iter->seq);
	}
}

void ContigData::copyUsableSeqSet(PSeqSet& outSeqs) const
{
	for(CKDataVec::const_iterator iter = m_kmerVec.begin(); iter != m_kmerVec.end(); ++iter)
	{
		if(iter->usable)
		{
			outSeqs.insert(iter->seq);
		}
	}
}
