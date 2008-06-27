#include "Aligner.h"

//
// Constructor
// 
Aligner::Aligner(int hashSize) : m_hashSize(hashSize)
{
	//m_pDatabase = new SeqPosHashMap(2 << 26);
	m_pDatabase = new SeqPosHashMap();
}

//
// Destructor
//
Aligner::~Aligner()
{
	delete m_pDatabase;
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
			p.contig = id;
			p.pos = i;
			m_pDatabase->insert(dbRecord(kmer, p));
		}
	}
}

void Aligner::alignRead(const PackedSeq& seq, AlignmentVector& alignVec)
{
	getAlignmentsInternal(seq, false, alignVec);
	getAlignmentsInternal(reverseComplement(seq), true, alignVec);
}

void Aligner::getAlignmentsInternal(const PackedSeq& seq, bool isRC, AlignmentVector& resultVector)
{
	// The results
	AlignmentSet aligns;

	int seqLen = seq.getSequenceLength();
	for(int i = 0; i < (seqLen - m_hashSize) + 1; ++i)
	{
		// Generate kmer
		PackedSeq kmer = seq.subseq(i, m_hashSize);

		// Get the alignment positions
		LookupResult result = m_pDatabase->equal_range(kmer);
		
		for (SPHMConstIter resultIter = result.first; resultIter != result.second; ++resultIter)
		{
			//printf("Seq: %s Contig: %s position: %d\n", seq.decode().c_str(), resultIter->second.contig.c_str(), resultIter->second.pos);
			int pos = resultIter->second.pos;
			ContigID contig = resultIter->second.contig;
			aligns[contig].insert(pos);
		}
	}
	
	coalesceAlignments(aligns, isRC, resultVector);
}

void Aligner::coalesceAlignments(const AlignmentSet& alignSet, bool isRC, AlignmentVector& resultVector)
{
	AlignmentResult result;
	AlignmentVector allAlignments;
	int bestLength = 0;
	
	// For each contig that this read hits, coalesce the alignments into contiguous groups
	for(AlignmentSet::const_iterator ctgIter = alignSet.begin(); ctgIter != alignSet.end(); ++ctgIter)
	{
		// Get the starting position
		assert(!ctgIter->second.empty());
		IntSet::iterator prevIter = ctgIter->second.begin();
		IntSet::iterator currIter = prevIter;
		currIter++;
		
		int currAlignStart = *prevIter;

		while(currIter != ctgIter->second.end())
		{
			// Discontinuity found
			if(*currIter != *prevIter + 1)
			{
				// create the alignment
				Alignment align = createAlignment(ctgIter->first, currAlignStart, *prevIter, isRC);
				
				// store the best length
				if(align.length > bestLength)
				{
					bestLength = align.length;
				}
				
				// save the alignment
				allAlignments.push_back(align);
					
				currAlignStart = *currIter;
			}
			
			prevIter = currIter;
			currIter++;
		}
		
		// Create the last alignment
		Alignment align = createAlignment(ctgIter->first, currAlignStart, *prevIter, isRC);
		if(align.length > bestLength)
		{
			bestLength = align.length;
		}		
		
		// save the alignment
		allAlignments.push_back(align);
	}
	
	// Filter the alignments to only return the best (longest) alignment(s)	
	// TODO: axe this? maybe return everything for super-strict alignments
	for(AlignmentVector::iterator alignIter = allAlignments.begin(); alignIter != allAlignments.end(); alignIter++)
	{
		//if(alignIter->length >= bestLength)
		{
			resultVector.push_back(*alignIter);	
		}	
	}
}

Alignment Aligner::createAlignment(ContigID contig, int start, int end, bool isRC)
{
	Alignment align;
	align.contig = contig;
	align.start = start;
	align.length = (end + m_hashSize) - start;
	align.isRC = isRC;
	return align;
}
