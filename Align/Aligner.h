#ifndef ALIGNER_H
#define ALIGNER_H 1

#include "config.h"
#include "Align/Options.h"
#include "ConstString.h"
#include "Functional.h"
#include "HashMap.h"
#include "Kmer.h"
#include <cassert>
#include <cstdlib>
#include <cstring> // for strcpy
#include <functional>
#include <iostream>
#include <istream>
#include <limits>
#include <map>
#include <ostream>
#include <string>
#include <vector>

typedef std::string StringID;

/** An ungapped alignment of a query to a target. */
struct Alignment
{
	StringID contig;
	int contig_start_pos;
	int read_start_pos;
	int align_length;
	int read_length;
	bool isRC;

	Alignment() { }

	Alignment(const Alignment& o,
			std::string /*qid*/, std::string /*seq*/) { *this = o; }

	Alignment(StringID contig, int contig_start, int read_start,
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

/** A tuple of a target ID and position. */
struct Position
{
	uint32_t contig;
	uint32_t pos; // 0 indexed
	Position(uint32_t contig = std::numeric_limits<uint32_t>::max(),
			uint32_t pos = std::numeric_limits<uint32_t>::max())
		: contig(contig), pos(pos) { }

	/** Mark this seed as a duplicate. */
	void setDuplicate(const char* thisContig, const char* otherContig,
			const Sequence& kmer)
	{
		if (opt::multimap == opt::IGNORE)
			contig = std::numeric_limits<uint32_t>::max();
		else {
			std::cerr << "error: duplicate k-mer in "
				<< thisContig
				<< " also in "
				<< otherContig
				<< ": " << kmer << '\n';
			exit(EXIT_FAILURE);
		}
	}

	/** Return whether this seed is a duplciate. */
	bool isDuplicate() const
	{
		return contig == std::numeric_limits<uint32_t>::max();
	}
};

typedef hash_multimap<Kmer, Position, hashKmer> SeqPosHashMultiMap;

#if HAVE_GOOGLE_SPARSE_HASH_MAP
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<Kmer, Position,
		hashKmer> SeqPosHashUniqueMap;
#else
typedef hash_map<Kmer, Position, hashKmer> SeqPosHashUniqueMap;
#endif


typedef std::vector<Alignment> AlignmentVector;

/**
 * Index a target sequence and align query sequences to that indexed
 * target.
 */
template <class SeqPosHashMap>
class Aligner
{
	public:
		typedef typename SeqPosHashMap::iterator map_iterator;
		typedef typename SeqPosHashMap::const_iterator
			map_const_iterator;

		Aligner(int hashSize, size_t buckets)
			: m_hashSize(hashSize), m_target(buckets) { }

		Aligner(int hashSize, size_t buckets, float factor)
			: m_hashSize(hashSize)
		{
			m_target.max_load_factor(factor);
			m_target.rehash(buckets);
		}

		void addReferenceSequence(const StringID& id,
				const Sequence& seq);
		void addReferenceSequence(const Kmer& kmer, Position pos);

		template <class oiterator>
		void alignRead(const std::string& qid, const Sequence& seq,
				oiterator dest);

		size_t size() const { return m_target.size(); }
		size_t bucket_count() const
		{
			return m_target.bucket_count();
		}

		/** Return the number of duplicate k-mer in the target. */
		size_t countDuplicates() const
		{
			assert(opt::multimap == opt::IGNORE);
			return count_if(m_target.begin(), m_target.end(),
					compose1(std::mem_fun_ref(&Position::isDuplicate),
						mem_var(&SeqPosHashMap::value_type::second)));
		}

	private:
		explicit Aligner(const Aligner&);

		typedef std::map<unsigned, AlignmentVector> AlignmentSet;

		void alignKmer(
				AlignmentSet& aligns, const Sequence& kmer,
				bool isRC, bool good, int read_ind, int seqLen);

		AlignmentSet getAlignmentsInternal(
				const Sequence& seq, bool isRC);

		template <class oiterator>
		void coalesceAlignments(
				const std::string& qid, const std::string& seq,
				const AlignmentSet& alignSet,
				oiterator& dest);

		// The number of bases to hash on
		int m_hashSize;

		/** A map of k-mer to contig coordinates. */
		SeqPosHashMap m_target;

		/** A dictionary of contig IDs. */
		std::vector<const_string> m_dict;

		unsigned contigIDToIndex(const std::string& id)
		{
			m_dict.push_back(id);
			return m_dict.size() - 1;
		}

		cstring contigIndexToID(unsigned index)
		{
			assert(index < m_dict.size());
			return m_dict[index];
		}
};

#endif
