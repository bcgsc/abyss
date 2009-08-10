#ifndef ALIGNER_H
#define ALIGNER_H 1

//
// Aligner - Simple class to do an approximate alignment of reads to the input reference sequence
//

#include "config.h"
#include "Dictionary.h"
#include "HashMap.h"
#include "PackedSeq.h"
#include <cassert>
#include <istream>
#include <map>
#include <ostream>
#include <string>
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

	Alignment() { }

	Alignment(ContigID contig, int contig_start, int read_start,
			int align_length, int read_length, bool isRC) :
		contig(contig),
		contig_start_pos(contig_start),
		read_start_pos(read_start),
		align_length(align_length),
		read_length(read_length),
		isRC(isRC)
	{
	}

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
		assert(tend <= tlength);
		rc.contig_start_pos = tlength - tend;
		rc.isRC = !isRC;
		return rc;
	}

	/** Return an alignment of the reverse complement of the query to
	 * the same target.
	 */
	Alignment flipQuery() const
	{
		Alignment rc(*this);
		unsigned qend = read_start_pos + align_length;
		assert(qend <= (unsigned)read_length);
		rc.read_start_pos = read_length - qend;
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

struct Position
{
	uint32_t contig;
	uint32_t pos; // 0 indexed
	Position(uint32_t contig, uint32_t pos)
		: contig(contig), pos(pos) { }
};

typedef hash_multimap<PackedSeq, Position,
		PackedSeqHasher, PackedSeqEqual> SeqPosHashMultiMap;

#if HAVE_GOOGLE_SPARSE_HASH_SET
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<PackedSeq, Position,
		PackedSeqHasher, PackedSeqEqual> SeqPosHashUniqueMap;
#else
typedef hash_map<PackedSeq, Position,
		PackedSeqHasher, PackedSeqEqual> SeqPosHashUniqueMap;
#endif


typedef std::vector<Alignment> AlignmentVector;
typedef std::map<unsigned, AlignmentVector> AlignmentSet;

template <class SeqPosHashMap>
class Aligner
{
	public:
		typedef typename SeqPosHashMap::const_iterator SPHMConstIter;
		typedef std::pair<SPHMConstIter, SPHMConstIter> LookupResult;

		Aligner(int hashSize, int buckets)
			: m_hashSize(hashSize), m_target(buckets) { }

		void addReferenceSequence(const ContigID& id, const Sequence& seq);

		// Align an individual sequence
		template <class oiterator>
		void alignRead(const Sequence& seq, oiterator dest);

		size_t size() const { return m_target.size(); }
		size_t bucket_count() const
		{
			return m_target.bucket_count();
		}

		/** Set the maximum load factor. */
		void max_load_factor(float factor)
		{
			m_target.max_load_factor(factor);
		}

	private:

		// Internal alignment function, perform the actual alignment
		template <class oiterator>
		void getAlignmentsInternal(const Sequence& seq, bool isRC,
				oiterator& dest);

		// Coalesce all the hash hits into contiguous alignments
		template <class oiterator>
		void coalesceAlignments(const AlignmentSet& alignSet,
				oiterator& dest);

		// The number of bases to hash on
		int m_hashSize;

		/** A map of k-mer to contig coordinates. */
		SeqPosHashMap m_target;

		/** A dictionary of contig IDs. */
		Dictionary m_contigDict;

		unsigned contigIDToIndex(const ContigID& id)
		{
			return m_contigDict.serial(id);
		}

		const ContigID& contigIndexToID(unsigned index)
		{
			return m_contigDict.key(index);
		}
};

#endif
