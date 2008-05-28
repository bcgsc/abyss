#ifndef ALIGNMENTCACHE_H
#define ALIGNMENTCACHE_H

#include "Aligner.h"
#include "PackedSeq.h"

//
// A class to hold kmer->contig lookup tables
//

typedef std::set<ContigID> ContigIDColl;
typedef std::map<PackedSeq, ContigIDColl> AlignDB;

class AlignmentCache
{
	public:
		AlignmentCache();
		
		void addAlignment(const PackedSeq& seq, const ContigID& id);
		
		//void serialize(std::string filename);
		//void unserialize(std::string filename);
		
	private:
		AlignDB m_alignmentCache;
	
};

#endif
