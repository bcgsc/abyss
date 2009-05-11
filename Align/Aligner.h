#ifndef ALIGNER_H
#define ALIGNER_H 1

//
// Aligner - Simple class to do an approximate alignment of reads to the input reference sequence
//

#include "config.h"
#include "PackedSeq.h"
#include <string>
#include <istream>
#include <ostream>
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
		unsigned tend = contig_start_pos + align_length;
		return !isRC ? contig_start_pos - read_start_pos
			: tend + read_start_pos;
	}

	/** Return the distance between the specified alignments.
	 * May be used to calculate fragment size when the alignments are
	 * mate pairs.
	 */
	int operator-(const Alignment& o) const
	{
		return targetAtQueryStart() - o.targetAtQueryStart();
	}

	/** This alignment is converted to the corresponding alignment of
	 * the same query to the reverse complement of the target.
	 */
	Alignment flipTarget(unsigned tlength) const
	{
		Alignment rc(*this);
		unsigned tend = contig_start_pos + align_length;
		rc.contig_start_pos = tlength - tend;
		rc.isRC = !isRC;
		return rc;
	}

	static int calculateReverseReadStart(int read_start_pos,
			int read_length, int align_length)
	{
		unsigned qend = read_start_pos + align_length;
		return read_length - qend;
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
	uint32_t pos; // 0 indexed
};

#if HAVE_GOOGLE_SPARSE_HASH_SET
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<PackedSeq,
		Position, PackedSeqHasher, PackedSeqEqual> SeqPosHashMap;
#else
# include "HashMap.h"
typedef hash_map<PackedSeq,
		Position, PackedSeqHasher, PackedSeqEqual> SeqPosHashMap;
#endif

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

		size_t size() const { return m_pDatabase->size(); }
		size_t bucket_count() const
		{
			return m_pDatabase->bucket_count();
		}

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
