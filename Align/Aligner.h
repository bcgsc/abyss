#ifndef ALIGNER_H
#define ALIGNER_H

//
// Aligner - Simple class to do an approximate alignment of reads to the input reference sequence
//

#include "CommonDefs.h"
#include "PackedSeq.h"
#include <ext/hash_map>

// Alignment information
struct Alignment
{
	ContigID contig;
	int start;
	int length;
	bool isRC;
};

// Typedef the database pairing
typedef std::pair<PackedSeq, Position> dbRecord;
typedef __gnu_cxx::hash_multimap<PackedSeq, Position, PackedSeqHasher, PackedSeqEqual> SeqPosHashMap;

typedef SeqPosHashMap::const_iterator SPHMConstIter;
typedef std::pair<SPHMConstIter, SPHMConstIter> LookupResult;

typedef std::set<int> IntSet;
typedef std::map<ContigID, IntSet > AlignmentSet;
typedef std::map<ContigID, int > AlignmentResult;

typedef std::vector<Alignment> AlignmentVector;



class Aligner
{
	public:
		Aligner(int hashSize);
		~Aligner();
		
		// Generate the database to align to
		void CreateDatabase(const ContigMap& refSeqs);
		
		// Align the vector of sequences
		void AlignReads(PSequenceVector seqs);
		
		// Align an individual sequence
		AlignmentVector GetAlignments(const PackedSeq& seq);		
		
	private:
	
		// Internal alignment function, perform the actual alignment
		void GetAlignmentsInternal(const PackedSeq& seq, bool isRC, AlignmentVector& resultVector);	
		
		// Coalesce all the hash hits into contiguous alignments
		void CoalesceAlignments(const AlignmentSet& alignSet, bool isRC, AlignmentVector& resultVector);
	
		// Create a single alignment from a start and end position
		Alignment CreateAlignment(ContigID contig, int start, int end, bool isRC);
		
		// The number of bases to hash on
		int m_hashSize;
		
		// The database of sequence hashes to alignment positions
		SeqPosHashMap* m_pDatabase;
	
};

#endif
