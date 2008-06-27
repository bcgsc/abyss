#ifndef ALIGNER_H
#define ALIGNER_H

//
// Aligner - Simple class to do an approximate alignment of reads to the input reference sequence
//

#include "CommonDefs.h"
#include "PackedSeq.h"
#include <ext/hash_map>
#include <sstream>
#include <iostream>

// Alignment information
struct Alignment
{
	ContigID contig;
	int start;
	int length;
	bool isRC;
	
	// input a record
	friend std::istream& operator>> (std::istream& in, Alignment& a)
	{
		// Read 1 record from the stream
		in >> a.contig;
		in >> a.start;
		in >> a.length;
		in >> a.isRC;
				
		return in;
	}
	
	// output a record
	friend std::ostream& operator<< (std::ostream& o, Alignment& a)
	{
		o << a.contig << " " << a.start << " " << a.length << " " << a.isRC;
		return o;
	}
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
		void addReferenceSequence(const ContigID& id, const Sequence& seq);
		
		// Align an individual sequence
		void alignRead(const PackedSeq& seq, AlignmentVector& alignVec);
		
		// Get the number of sequences in the database
		size_t getNumSeqs() const { return m_pDatabase->size(); }
		
	private:
	
		// Internal alignment function, perform the actual alignment
		void getAlignmentsInternal(const PackedSeq& seq, bool isRC, AlignmentVector& resultVector);	
		
		// Coalesce all the hash hits into contiguous alignments
		void coalesceAlignments(const AlignmentSet& alignSet, bool isRC, AlignmentVector& resultVector);
	
		// Create a single alignment from a start and end position
		Alignment createAlignment(ContigID contig, int start, int end, bool isRC);
		
		// The number of bases to hash on
		int m_hashSize;
		
		// The database of sequence hashes to alignment positions
		SeqPosHashMap* m_pDatabase;
	
};

#endif
