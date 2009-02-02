#ifndef ALIGNER_H
#define ALIGNER_H

//
// Aligner - Simple class to do an approximate alignment of reads to the input reference sequence
//

#include "PackedSeq.h"
#include <ext/hash_map>
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

typedef std::string ContigID;

// Alignment information
struct Alignment
{
	ContigID contig;
	int contig_start_pos;
	int read_start_pos;
	int align_length;
	int read_length;
	bool isRC;

	/**
	 * Return the taret position at the query start.
	 * Note: not alignment start, and may be negative
	 */
	int targetAtQueryStart() const
	{
		unsigned e = read_start_pos + align_length;
		unsigned s = !isRC ? read_start_pos : read_length - e;
		return contig_start_pos - s;
	}

	/**
	 * Return the target position at the query end.
	 * Note: not alignment end
	 */
	int targetAtQueryEnd() const
	{
		return targetAtQueryStart() + read_length;
	}

	// flip the alignment with respect to the contig size
	void flip(int contigLength)
	{
		// flip the contig start position
		contig_start_pos = contigLength - (contig_start_pos + align_length);
		
		// flip the read start pos
		read_start_pos = calculateReverseReadStart(read_start_pos, read_length, align_length);
		
		// flip the rc bit
		isRC = !isRC;
	}
	
	static int calculateReverseReadStart(int pos, int readLen, int alignLength)
	{
		return (readLen - pos) - alignLength;
	}
	
	// input a record
	friend std::istream& operator>> (std::istream& in, Alignment& a)
	{
		// Read 1 record from the stream
		in >> a.contig;
		in >> a.contig_start_pos;
		in >> a.read_start_pos;
		in >> a.align_length;
		in >> a.read_length;
		in >> a.isRC;
				
		return in;
	}
	
	// output a record
	friend std::ostream& operator<< (std::ostream& o, const Alignment& a)
	{
		o << a.contig << " " << a.contig_start_pos << " " << a.read_start_pos << " " << a.align_length << " " << a.read_length << " " << a.isRC;
		return o;
	}
};

static inline int compareContigPos(const Alignment& a1, const Alignment& a2)
{
	return a1.contig_start_pos < a2.contig_start_pos;
}


struct Position
{
	ContigID contig;
	int pos; // 0 indexed
};

// Typedef the database pairing
typedef std::pair<PackedSeq, Position> dbRecord;
typedef __gnu_cxx::hash_multimap<PackedSeq, Position, PackedSeqHasher, PackedSeqEqual> SeqPosHashMap;

typedef SeqPosHashMap::const_iterator SPHMConstIter;
typedef std::pair<SPHMConstIter, SPHMConstIter> LookupResult;

typedef std::vector<Alignment> AlignmentVector;
typedef std::map<ContigID, AlignmentVector > AlignmentSet;
typedef std::map<ContigID, int > AlignmentResult;





class Aligner
{
	public:
		Aligner(int hashSize);
		~Aligner();
		
		// Generate the database to align to
		void addReferenceSequence(const ContigID& id, const Sequence& seq);
		
		// Align an individual sequence
		void alignRead(const Sequence& seq, AlignmentVector& alignVec);
		
		// Get the number of sequences in the database
		size_t getNumSeqs() const { return m_pDatabase->size(); }
		
	private:
	
		// Internal alignment function, perform the actual alignment
		void getAlignmentsInternal(const Sequence& seq, bool isRC, AlignmentVector& resultVector);	
		
		// Coalesce all the hash hits into contiguous alignments
		void coalesceAlignments(const AlignmentSet& alignSet, bool isRC, AlignmentVector& resultVector);
	
		// Create a single alignment from a start and end position
		Alignment createAlignment(ContigID contig, int contig_start, int read_start, int align_length, int read_length, bool isRC);
		
		// The number of bases to hash on
		int m_hashSize;
		
		// The database of sequence hashes to alignment positions
		SeqPosHashMap* m_pDatabase;
	
};

#endif
