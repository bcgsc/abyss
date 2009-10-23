#include "Aligner.h"
#include "Align/Options.h"
#include "PrefixIterator.h"
#include "Sequence.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <utility>

using namespace std;

namespace opt {
	/** For a duplicate k-mer in the target
	 * ERROR: report an error and exit
	 * MULTIMAP: report all alignments
	 * IGNORE: do not report any alignments
	 */
	int multimap;
};

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
	pair<map_iterator, bool> inserted
		= m_target.insert(make_pair(kmer, pos));
	if (inserted.second)
		return;
	if (opt::multimap == opt::IGNORE) {
		inserted.first->second.setDuplicate();
	} else {
		cerr << "error: duplicate k-mer in "
			<< contigIndexToID(inserted.first->second.contig)
			<< " also in "
			<< contigIndexToID(pos.contig)
			<< ": " << kmer.decode() << '\n';
		exit(EXIT_FAILURE);
	}
}

/** Create an index of the target sequence. */
template <class SeqPosHashMap>
void Aligner<SeqPosHashMap>::addReferenceSequence(const ContigID& id, const Sequence& seq)
{
	// Break the ref sequence into kmers of the hash size
	int size = seq.length();
	for(int i = 0; i < (size - m_hashSize + 1); ++i)
	{
		Sequence subseq(seq, i, m_hashSize);
		if (subseq.find("N") != string::npos)
			continue;
		addReferenceSequence(Kmer(subseq),
				Position(contigIDToIndex(id), i));
	}
}

template <class SeqPosHashMap>
template <class oiterator>
void Aligner<SeqPosHashMap>::alignRead(const Sequence& seq,
		oiterator dest)
{
	getAlignmentsInternal(seq, false, dest);
	getAlignmentsInternal(reverseComplement(seq), true, dest);
}

template <class SeqPosHashMap>
template <class oiterator>
void Aligner<SeqPosHashMap>::
getAlignmentsInternal(const Sequence& seq, bool isRC,
		oiterator& dest)
{
	// The results
	AlignmentSet aligns;

	int seqLen = seq.length();
	for(int i = 0; i < (seqLen - m_hashSize) + 1; ++i)
	{
		pair<map_const_iterator, map_const_iterator> range
			= m_target.equal_range(Kmer(seq.substr(i, m_hashSize)));
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

	coalesceAlignments(aligns, dest);
}

static int compareQueryPos(const Alignment& a1, const Alignment& a2)
{
	return a1.read_start_pos < a2.read_start_pos;
}

/** Coalesce the k-mer alignments into a read alignment. */
template <class SeqPosHashMap>
template <class oiterator>
void Aligner<SeqPosHashMap>::
coalesceAlignments(const AlignmentSet& alignSet, oiterator& dest)
{
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
				*dest++ = currAlign;
				currAlign = *currIter;
			} else {
				currAlign.align_length++;
				if (currAlign.isRC)
					currAlign.contig_start_pos--;
			}

			prevIter = currIter;
			currIter++;
		}

		*dest++ = currAlign;
	}
}

// Explicit instantiation.
template void Aligner<SeqPosHashMultiMap>::addReferenceSequence(
		const ContigID& id, const Sequence& seq);

template void Aligner<SeqPosHashUniqueMap>::addReferenceSequence(
		const ContigID& id, const Sequence& seq);

template void Aligner<SeqPosHashMultiMap>::
alignRead<prefix_ostream_iterator<Alignment> >(
		const Sequence& seq, prefix_ostream_iterator<Alignment> dest);

template void Aligner<SeqPosHashUniqueMap>::
alignRead<prefix_ostream_iterator<Alignment> >(
		const Sequence& seq, prefix_ostream_iterator<Alignment> dest);
