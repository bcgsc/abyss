#ifndef ALIGNMENTCACHE_H
#define ALIGNMENTCACHE_H

#include "Aligner.h"
#include "PackedSeq.h"

//
// A class to hold kmer->contig lookup tables
//

typedef std::set<ContigID> ContigIDSet;
typedef std::map<PackedSeq, ContigIDSet> AlignDB;

class AlignmentCache
{
	public:
		AlignmentCache();
		
		void addKeys(const PSeqSet& seqSet, const ContigID& id);
		void deleteKeys(const PSeqSet& seqSet, const ContigID& id);
		void getSet(const PackedSeq& seq, ContigIDSet& outset) const;
		
		bool compare(const AlignmentCache& otherDB);
		
		void printSet(const ContigIDSet& seqSet) const;
		
		void concatSets(ContigIDSet& seqSet1, const ContigIDSet& seqSet2) const;
		
		//void serialize(std::string filename);
		//void unserialize(std::string filename);
		
	private:
		void addAlignment(const PackedSeq& seq, const ContigID& id);
		void removeAlignment(const PackedSeq& seq, const ContigID& id);		
		AlignDB m_alignmentCache;
	
};

#endif
