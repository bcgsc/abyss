#ifndef CONTIGDATA_H
#define CONTIGDATA_H

#include "CommonDefs.h"
#include "PackedSeq.h"

struct ContigKmerData
{
	PackedSeq seq;
	bool usable;
};

typedef std::vector<ContigKmerData> CKDataVec;

class ContigData
{
	public:
		
		ContigData(const ContigID& id, const Sequence& s, size_t kmer);
		
		void addID(const ContigID& id) { m_idParts.insert(id); }
		
		void addSeqs(const Sequence& s, CKDataVec::iterator position, bool usable);
		
		void copySeqSet(PSeqSet& outSeqs) const;
		void copyUsableSeqSet(PSeqSet& outSeqs) const;
		
		ContigIDSet getExclusionSet() const { return m_idParts; }
		
		void validateSequences() const;
		
		Sequence m_seq;
		CKDataVec m_kmerVec;
		
		// The record of initial contigs that made up this contig
		ContigIDSet m_idParts;
		
		size_t m_kmer;

	
	
};

#endif
