#ifndef ALIGNER_H
#define ALIGNER_H 1

//
// Aligner - Simple class to do an approximate alignment of reads to the input reference sequence
//

#include "config.h"
#include "Dictionary.h"
#include "HashMap.h"
#include "Kmer.h"
#include <cassert>
#include <istream>
#include <map>
#include <ostream>
#include <string>
#include <vector>
#include <limits> //for uint32_t.max

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

	friend std::istream& operator >>(std::istream& in, Alignment& a)
	{
		return in >> a.contig
			>> a.contig_start_pos
			>> a.read_start_pos
			>> a.align_length
			>> a.read_length
			>> a.isRC;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const Alignment& a)
	{
		return out << a.contig << ' '
			<< a.contig_start_pos << ' '
			<< a.read_start_pos << ' '
			<< a.align_length << ' '
			<< a.read_length << ' '
			<< a.isRC;
	}
};

struct Position
{
	uint32_t contig;
	uint32_t pos; // 0 indexed
	Position(uint32_t contig, uint32_t pos)
		: contig(contig), pos(pos) { }
	void setDuplicate() {contig = std::numeric_limits<uint32_t>::max ();}
	bool isDuplicate() const {return contig == std::numeric_limits<uint32_t>::max();}
};

struct KmerEqual
{
	bool operator()(const Kmer& a, const Kmer& b) const
	{
		return a == b;
	}
};

struct KmerHasher
{
	size_t operator()(const Kmer& o) const { return o.getHashCode(); }
};

typedef hash_multimap<Kmer, Position,
		KmerHasher, KmerEqual> SeqPosHashMultiMap;

#if HAVE_GOOGLE_SPARSE_HASH_SET
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<Kmer, Position,
		KmerHasher, KmerEqual> SeqPosHashUniqueMap;
#else
typedef hash_map<KmerSeq, Position,
		KmerSeqHasher, KmerSeqEqual> SeqPosHashUniqueMap;
#endif


typedef std::vector<Alignment> AlignmentVector;
typedef std::map<unsigned, AlignmentVector> AlignmentSet;

template <class SeqPosHashMap>
class Aligner
{
	public:
		typedef typename SeqPosHashMap::iterator map_iterator;
		typedef typename SeqPosHashMap::const_iterator
			map_const_iterator;

		Aligner(int hashSize, int buckets)
			: m_hashSize(hashSize), m_target(buckets) { }

		Aligner(int hashSize, int buckets, float factor)
			: m_hashSize(hashSize)
		{
			m_target.max_load_factor(factor);
			m_target.rehash(buckets);
		}

		void addReferenceSequence(const ContigID& id,
				const Sequence& seq);
		void addReferenceSequence(const Kmer& kmer, Position pos);

		template <class oiterator>
		void alignRead(const Sequence& seq, oiterator dest);

		size_t size() const { return m_target.size(); }
		size_t bucket_count() const
		{
			return m_target.bucket_count();
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
