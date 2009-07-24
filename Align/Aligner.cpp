#include "Aligner.h"
#include "PrefixIterator.h"
#include "Sequence.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <utility>

using namespace std;

//
// Constructor
//
Aligner::Aligner(int hashSize) : m_hashSize(hashSize)
{
#if HAVE_GOOGLE_SPARSE_HASH_SET
	m_pDatabase = new SeqPosHashMap(1 << 28);
	m_pDatabase->max_load_factor(0.2);
#else
	m_pDatabase = new SeqPosHashMap(1 << 26);
#endif
}

//
// Destructor
//
Aligner::~Aligner()
{
	delete m_pDatabase;
}

/** Convert the specified contig ID string to a numeric index. */
unsigned Aligner::contigIDToIndex(ContigID id)
{
	pair<ContigDict::const_iterator, bool> inserted
		= m_contigDict.insert(make_pair(id, m_contigDict.size()));
	if (inserted.second)
		m_contigIDs.push_back(id);
	return inserted.first->second;
}

/** Convert the specified contig ID numeric index to a string. */
const ContigID& Aligner::contigIndexToID(unsigned index)
{
	return m_contigIDs.at(index);
}

//
// Create the database for the reference sequences
//
void Aligner::addReferenceSequence(const ContigID& id, const Sequence& seq)
{		
	// Break the ref sequence into kmers of the hash size
	int size = seq.length();
	for(int i = 0; i < (size - m_hashSize + 1); ++i)
	{
		Sequence subseq = seq.substr(i, m_hashSize);
		
		// skip seqs that have unknown bases
		if(subseq.find("N") == std::string::npos)
		{
			PackedSeq kmer(subseq);
			//printf("indexed seq: %s\n", kmer.decode().c_str());
			Position p;
			p.contig = contigIDToIndex(id);
			p.pos = i;
#if HAVE_GOOGLE_SPARSE_HASH_SET
			if (m_pDatabase->count(kmer) > 0) {
				cerr << "error: duplicate k-mer in "
					<< (*m_pDatabase)[kmer].contig << " also in "
					<< id << ": " << kmer.decode() << '\n';
				exit(EXIT_FAILURE);
			}
#endif
			m_pDatabase->insert(make_pair(kmer, p));
		}
	}
}

template<typename oiterator>
void Aligner::alignRead(const Sequence& seq,
		oiterator dest)
{
	getAlignmentsInternal(seq, false, dest);
	getAlignmentsInternal(reverseComplement(seq), true, dest);
}

// Explicit instantiation.
template void Aligner::alignRead
		<prefix_ostream_iterator<Alignment> >(
		const Sequence& seq,
		prefix_ostream_iterator<Alignment> dest);

template<typename oiterator>
void Aligner::getAlignmentsInternal(const Sequence& seq, bool isRC,
		oiterator& dest)
{
	// The results
	AlignmentSet aligns;

	int seqLen = seq.length();
	for(int i = 0; i < (seqLen - m_hashSize) + 1; ++i)
	{
		// Generate kmer
		PackedSeq kmer = seq.substr(i, m_hashSize);

		// Get the alignment positions
		LookupResult result = m_pDatabase->equal_range(kmer);
		
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

template<typename oiterator>
void Aligner::coalesceAlignments(const AlignmentSet& alignSet,
		oiterator& dest)
{
	// For each contig that this read hits, coalesce the alignments into contiguous groups
	for(AlignmentSet::const_iterator ctgIter = alignSet.begin(); ctgIter != alignSet.end(); ++ctgIter)
	{
		AlignmentVector alignVec = ctgIter->second;
		
		if(alignVec.empty())
		{
			continue;
		}
		
		// Sort the alignment vector by contig alignment position
		sort(alignVec.begin(), alignVec.end(), compareContigPos);
				
		// Get the starting position
		assert(!ctgIter->second.empty());
		AlignmentVector::iterator prevIter = alignVec.begin();
		AlignmentVector::iterator currIter = alignVec.begin() + 1;
		
		Alignment currAlign = *prevIter;

		while(currIter != alignVec.end())
		{
			//std::cout << "CurrAlign: " << currAlign << "\n";
			//std::cout << "AlignIter: " << *currIter << "\n";
			
			// Discontinuity found
			if(currIter->contig_start_pos != prevIter->contig_start_pos + 1)
			{
				//std::cout << "	Discontinous, saving\n";
				*dest++ = currAlign;
				currAlign = *currIter;
			}
			else
			{
				//bstd::cout << "	Continous, updating\n";
				// alignments are consistent, increase the length of the alignment
				currAlign.align_length++;
				currAlign.read_start_pos = std::min(currAlign.read_start_pos, currIter->read_start_pos);
			}
			
			prevIter = currIter;
			currIter++;
		}

		*dest++ = currAlign;
	}
}
