#include "Aligner.h"
#include "Iterator.h"
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
					contigIndexToID(pos.contig), kmer.str());
		return;
	}

	pair<map_iterator, bool> inserted
		= m_target.insert(make_pair(kmer, pos));

	if (!inserted.second)
		if (!inserted.first->second.isDuplicate())
			inserted.first->second.setDuplicate(
					contigIndexToID(inserted.first->second.contig),
					contigIndexToID(pos.contig), kmer.str());
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

/** Store all alignments for a given Kmer in the parameter aligns.
 *  @param[out] aligns Map of contig IDs to alignment vectors.
 */
template <class SeqPosHashMap>
void Aligner<SeqPosHashMap>::alignKmer(
		AlignmentSet& aligns, const Sequence& seq,
		bool isRC, bool good, int read_ind, int seqLen)
{
	assert(read_ind >= 0);
	Sequence kmer(seq, read_ind, m_hashSize);
	if (!good && kmer.find_first_not_of("ACGT0123") != string::npos)
		return;

	pair<map_const_iterator, map_const_iterator> range
		= m_target.equal_range(Kmer(kmer));

	if (range.first != range.second
				&& opt::multimap == opt::IGNORE
				&& range.first->second.isDuplicate())
		return;

	for (map_const_iterator resultIter = range.first;
			resultIter != range.second; ++resultIter) {
		assert(opt::multimap != opt::IGNORE
				|| !resultIter->second.isDuplicate());

		int read_pos = !isRC ? read_ind
			: Alignment::calculateReverseReadStart(
					read_ind, seqLen, m_hashSize);
		unsigned ctgIndex = resultIter->second.contig;
		Alignment align(string(),
				resultIter->second.pos, read_pos, m_hashSize,
				seqLen, isRC);
		aligns[ctgIndex].push_back(align);
	}
}

template <class SeqPosHashMap>
typename Aligner<SeqPosHashMap>::AlignmentSet
Aligner<SeqPosHashMap>::getAlignmentsInternal(
		const Sequence& seq, bool isRC)
{
	// The results
	AlignmentSet aligns;

	bool good = seq.find_first_not_of("ACGT0123") == string::npos;
	int seqLen = seq.length();
	int last_kmer = seqLen - m_hashSize;

	if (last_kmer < 0)
		return aligns;

	// Align the first kmer
	alignKmer(aligns, seq, isRC, good, 0, seqLen);

	if (last_kmer == 0)
		return aligns;

	// Align the last kmer
	alignKmer(aligns, seq, isRC, good, last_kmer, seqLen);

	// Short-cut logic ignoring the middle alignments if the first
	// and last kmers overlap, and align to the same contig
	if (good && seqLen <= 2 * m_hashSize && aligns.size() == 1) {
		AlignmentSet::const_iterator ctgIter = aligns.begin();
		const AlignmentVector& a = ctgIter->second;
		if (ctgIter->second.size() == 2) {
			int qstep = isRC
					? a[0].read_start_pos - a[1].read_start_pos
					: a[1].read_start_pos - a[0].read_start_pos;
			assert(qstep >= 0);

			// Verify this isn't a kmer aligning to two parts of a
			// contig, and that the alignments are coalescable.
			if (qstep == last_kmer &&
					a[1].contig_start_pos
						== a[0].contig_start_pos + qstep)
				return aligns;
		}
	}

	// Align middle kmers
	for(int i = 1; i < last_kmer; ++i)
		alignKmer(aligns, seq, isRC, good, i, seqLen);

	return aligns;
}

static int compareQueryPos(const Alignment& a1, const Alignment& a2)
{
	return a1.read_start_pos < a2.read_start_pos;
}

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
			int qstep = currIter->read_start_pos -
						prevIter->read_start_pos;
			assert(qstep >= 0 && qstep < m_hashSize);
			int tstep = currIter->isRC ? -qstep : qstep;
			if (currIter->contig_start_pos
						== prevIter->contig_start_pos + tstep
					&& qstep < m_hashSize) {
				currAlign.align_length += qstep;
				if (currAlign.isRC)
					currAlign.contig_start_pos -= qstep;
			} else {
				currAlign.contig = contigIndexToID(ctgIter->first);
				*dest++ = value_type(currAlign, qid, seq);
				currAlign = *currIter;
			}

			prevIter = currIter;
			currIter++;
		}

		currAlign.contig = contigIndexToID(ctgIter->first);
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
