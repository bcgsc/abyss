#ifndef ALIGNMENTCACHE_H
#define ALIGNMENTCACHE_H

#include "Aligner.h"
#include "PackedSeq.h"

struct AlignData
{
	ContigID contigID;
	int position;
	bool isRC;
	
	bool operator<(const AlignData& a2) const
	{
		int cmp = strcmp(this->contigID.c_str(), a2.contigID.c_str());
		if(cmp == 0)
		{
			// strings are equal, compare on position
			return (this->position < a2.position);
		}
		else
		{
			return cmp < 0;
		}
	}
	
	friend std::ostream& operator<<(std::ostream& out, const AlignData& object)
	{
		out << "(" << object.contigID << "," << object.position << "," << object.isRC << ")";
		return out;
	}	
};

//
// A class to hold kmer->contig lookup tables
//

typedef std::set<AlignData> AlignSet;
typedef std::map<PackedSeq, AlignSet> AlignDB;

class AlignmentCache
{
	public:
		AlignmentCache();
		
		void addAlignment(const PackedSeq& seq, const ContigID& id, int position);
		void removeAlignment(const PackedSeq& seq, const ContigID& id, int position);
		void getAlignments(const PackedSeq& seq, AlignSet& outset) const;
		
		void printAlignmentsForSeq(const PackedSeq& seq) const;
		void translatePosition(const PackedSeq& seq, const ContigID& id, int oldPosition, int newPosition);
		
		/*
		void addKeys(const PSeqSet& seqSet, const ContigID& id);
		void deleteKeys(const PSeqSet& seqSet, const ContigID& id);
		*/
		void getSet(const PackedSeq& seq, ContigIDSet& outset) const;
		
		bool compare(const AlignmentCache& otherDB);

		void concatSets(ContigIDSet& seqSet1, const ContigIDSet& seqSet2) const;
	
		//void serialize(std::string filename);
		//void unserialize(std::string filename);
		
	private:

		void removeAlignmentInternal(const PackedSeq& seq, const ContigID& id, int position);
		AlignDB m_alignmentCache;
	
};

#endif
