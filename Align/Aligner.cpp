#include "Aligner.h"
#include "PrefixIterator.h"
#include "Sequence.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <utility>

using namespace std;

namespace opt {
	extern int multimap;
};

/** Create an index of the target sequence. */
template <class SeqPosHashMap>
void Aligner<SeqPosHashMap>::addReferenceSequence(const ContigID& id, const Sequence& seq)
{
	// Break the ref sequence into kmers of the hash size
	int size = seq.length();
	for(int i = 0; i < (size - m_hashSize + 1); ++i)
	{
		Sequence subseq = seq.substr(i, m_hashSize);
		if (subseq.find("N") != string::npos)
			continue;

		PackedSeq kmer(subseq);
		if (!opt::multimap) {
			class SeqPosHashMap::const_iterator it
				= m_target.find(kmer);
			if (it != m_target.end()) {
				cerr << "error: duplicate k-mer in "
					<< contigIndexToID(it->second.contig)
					<< " also in " << id << ": "
					<< kmer.decode() << '\n';
				exit(EXIT_FAILURE);
			}
		}
		m_target.insert(make_pair(kmer,
					Position(contigIDToIndex(id), i)));
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
		PackedSeq kmer = seq.substr(i, m_hashSize);
		LookupResult result = m_target.equal_range(kmer);

		for (SPHMConstIter resultIter = result.first; resultIter != result.second; ++resultIter)
		{
			//printf("Seq: %s Contig: %s position: %d\n", seq.decode().c_str(), resultIter->second.contig.c_str(), resultIter->second.pos);
			int read_pos;
			
			// The read position coordinate is wrt to the forward read position
			if(!isRC)
			{
				read_pos = i;
			}
			else
			{
				read_pos = Alignment::calculateReverseReadStart(i, seqLen, m_hashSize);
			}

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
