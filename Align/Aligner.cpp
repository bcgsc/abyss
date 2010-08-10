#include "Aligner.h"
#include "AffixIterator.h"
#include "SAM.h"
#include "Sequence.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iterator>
#include <utility>

using namespace std;

namespace opt {
	/** For a duplicate k-mer in the target
	 * ERROR: report an error and exit
	 * MULTIMAP: report all alignments
	 * IGNORE: do not report any alignments
	 */
	int multimap;
}

template <>
void Aligner<SeqPosHashMultiMap>::addReferenceSequence(
		const Kmer& kmer, Position pos)
{
	assert(opt::multimap == opt::MULTIMAP);
	m_target.insert(make_pair(kmer, pos));
}

template <class SeqPosHashMap>
void Aligner<SeqPosHashMap>::addReferenceSequence(
		const Kmer& kmer, Position pos)
{
	assert(opt::multimap != opt::MULTIMAP);
	Kmer rc_kmer = reverseComplement(kmer);
	map_iterator findIt = m_target.find(rc_kmer);

	if (findIt != m_target.end()) {
		if (!findIt->second.isDuplicate())
			findIt->second.setDuplicate(
					contigIndexToID(findIt->second.contig),
					contigIndexToID(pos.contig), kmer.decode());
		return;
	}

	pair<map_iterator, bool> inserted
		= m_target.insert(make_pair(kmer, pos));

	if (!inserted.second)
		if (!inserted.first->second.isDuplicate())
			inserted.first->second.setDuplicate(
					contigIndexToID(inserted.first->second.contig),
					contigIndexToID(pos.contig), kmer.decode());
}

/** Create an index of the target sequence. */
template <class SeqPosHashMap>
void Aligner<SeqPosHashMap>::addReferenceSequence(
		const StringID& idString, const Sequence& seq)
{
	unsigned id = contigIDToIndex(idString);
	int size = seq.length();
	for(int i = 0; i < (size - m_hashSize + 1); ++i)
	{
		Sequence subseq(seq, i, m_hashSize);
		if (subseq.find("N") != string::npos)
			continue;
		addReferenceSequence(Kmer(subseq), Position(id, i));
	}
}

template <class SeqPosHashMap>
template <class oiterator>
void Aligner<SeqPosHashMap>::alignRead(
		const string& qid, const Sequence& seq,
		oiterator dest)
{
	coalesceAlignments(qid, seq,
			getAlignmentsInternal(seq, false), dest);
	Sequence seqrc = reverseComplement(seq);
	coalesceAlignments(qid, seqrc,
			getAlignmentsInternal(seqrc, true), dest);
}

template <class SeqPosHashMap>
typename Aligner<SeqPosHashMap>::AlignmentSet
Aligner<SeqPosHashMap>::getAlignmentsInternal(
		const Sequence& seq, bool isRC)
{
	// The results
	AlignmentSet aligns;

	bool good = seq.find_first_not_of("ACGT0123") == string::npos;
	bool discarded = true;
	int seqLen = seq.length();
	for(int i = 0; i < (seqLen - m_hashSize) + 1; ++i)
	{
		Sequence kmer(seq, i, m_hashSize);
		if (!good && kmer.find_first_not_of("ACGT0123")
				!= string::npos)
			continue;
		discarded = false;
		pair<map_const_iterator, map_const_iterator> range
			= m_target.equal_range(Kmer(kmer));
		for (map_const_iterator resultIter = range.first;
				resultIter != range.second; ++resultIter) {
			if (opt::multimap == opt::IGNORE
					&& resultIter->second.isDuplicate())
				break;

			int read_pos = !isRC ? i
				: Alignment::calculateReverseReadStart(
						i, seqLen, m_hashSize);
			unsigned ctgIndex = resultIter->second.contig;
			Alignment align(contigIndexToID(ctgIndex),
					resultIter->second.pos, read_pos, m_hashSize,
					seqLen, isRC);
			aligns[ctgIndex].push_back(align);
		}
	}
	return aligns;
}

static int compareQueryPos(const Alignment& a1, const Alignment& a2)
{
	return a1.read_start_pos < a2.read_start_pos;
}

/** Traits of an output iterator. */
template<class OutIter>
struct output_iterator_traits {
	typedef typename OutIter::value_type value_type;
};

/** Traits of a ostream_iterator. */
template<class T, class charT, class traits>
struct output_iterator_traits
	<std::ostream_iterator<T, charT, traits> >
{
	typedef T value_type;
};

/** Coalesce the k-mer alignments into a read alignment. */
template <class SeqPosHashMap>
template <class oiterator>
void Aligner<SeqPosHashMap>::coalesceAlignments(
		const string& qid, const string& seq,
		const AlignmentSet& alignSet,
		oiterator& dest)
{
	typedef typename output_iterator_traits<oiterator>::value_type
		value_type;
	for (AlignmentSet::const_iterator ctgIter = alignSet.begin();
			ctgIter != alignSet.end(); ++ctgIter) {
		AlignmentVector alignVec = ctgIter->second;
		assert(!alignVec.empty());

		sort(alignVec.begin(), alignVec.end(), compareQueryPos);

		AlignmentVector::iterator prevIter = alignVec.begin();
		AlignmentVector::iterator currIter = alignVec.begin() + 1;
		Alignment currAlign = *prevIter;
		while (currIter != alignVec.end()) {
			int qstep = 1;
			int tstep = currIter->isRC ? -1 : 1;
			if (currIter->read_start_pos
						!= prevIter->read_start_pos + qstep
					|| currIter->contig_start_pos
						!= prevIter->contig_start_pos + tstep) {
				*dest++ = value_type(currAlign, qid, seq);
				currAlign = *currIter;
			} else {
				currAlign.align_length++;
				if (currAlign.isRC)
					currAlign.contig_start_pos--;
			}

			prevIter = currIter;
			currIter++;
		}

		*dest++ = value_type(currAlign, qid, seq);
	}
}

// Explicit instantiation.
template void Aligner<SeqPosHashMultiMap>::addReferenceSequence(
		const StringID& id, const Sequence& seq);

template void Aligner<SeqPosHashUniqueMap>::addReferenceSequence(
		const StringID& id, const Sequence& seq);

template void Aligner<SeqPosHashMultiMap>::
alignRead<affix_ostream_iterator<Alignment> >(
		const string& qid, const Sequence& seq,
		affix_ostream_iterator<Alignment> dest);

template void Aligner<SeqPosHashUniqueMap>::
alignRead<affix_ostream_iterator<Alignment> >(
		const string& qid, const Sequence& seq,
		affix_ostream_iterator<Alignment> dest);

template void Aligner<SeqPosHashMultiMap>::
alignRead<ostream_iterator<SAMRecord> >(
		const string& qid, const Sequence& seq,
		ostream_iterator<SAMRecord> dest);

template void Aligner<SeqPosHashUniqueMap>::
alignRead<ostream_iterator<SAMRecord> >(
		const string& qid, const Sequence& seq,
		ostream_iterator<SAMRecord> dest);
