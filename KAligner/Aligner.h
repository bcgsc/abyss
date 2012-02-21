#ifndef ALIGNER_H
#define ALIGNER_H 1

#include "config.h"
#include "KAligner/Options.h"
#include "Alignment.h"
#include "ConstString.h"
#include "Functional.h"
#include "Kmer.h"
#include "UnorderedMap.h"
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

typedef unordered_multimap<Kmer, Position, hashKmer>
	SeqPosHashMultiMap;

#if HAVE_GOOGLE_SPARSE_HASH_MAP
# include <google/sparse_hash_map>
typedef google::sparse_hash_map<Kmer, Position,
		hashKmer> SeqPosHashUniqueMap;
#else
typedef unordered_map<Kmer, Position, hashKmer> SeqPosHashUniqueMap;
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
