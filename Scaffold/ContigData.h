#ifndef CONTIGDATA_H
#define CONTIGDATA_H

#include "CommonDefs.h"
#include "PackedSeq.h"

class ContigData
{
	public:
		ContigData(const Sequence& s, size_t kmer);
		
		Sequence seq;
		PSeqSet m_seqSet;
		size_t m_kmer;
	
	
};

#endif
