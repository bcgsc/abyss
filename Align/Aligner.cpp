#include "Aligner.h"
#include "Sequence.h"

//
// Constructor
// 
Aligner::Aligner(int hashSize) : m_hashSize(hashSize)
{
	//m_pDatabase = new SeqPosHashMap(2 << 26);
	m_pDatabase = new SeqPosHashMap(2 << 25);
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

void Aligner::alignRead(const Sequence& seq, AlignmentVector& alignVec)
{
	getAlignmentsInternal(seq, false, alignVec);
	getAlignmentsInternal(reverseComplement(seq), true, alignVec);
}

void Aligner::getAlignmentsInternal(const Sequence& seq, bool isRC, AlignmentVector& resultVector)
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
			
			Alignment align = createAlignment(resultIter->second.contig, resultIter->second.pos, read_pos, m_hashSize, seqLen, isRC);
			aligns[resultIter->second.contig].push_back(align);
		}
	}
	
	coalesceAlignments(aligns, isRC, resultVector);
}

void Aligner::coalesceAlignments(const AlignmentSet& alignSet, bool /*isRC*/, AlignmentVector& resultVector)
{
	AlignmentResult result;
	AlignmentVector allAlignments;
	
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
				// save the alignment
				resultVector.push_back(currAlign);	
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
		
		// save the alignment
		resultVector.push_back(currAlign);
	}
	
}

Alignment Aligner::createAlignment(ContigID contig, int contig_start, int read_start, int align_length, int read_length, bool isRC)
{
	Alignment align;
	align.contig = contig;
	align.contig_start_pos = contig_start;
	align.read_start_pos = read_start;
	align.align_length = align_length;
	align.read_length = read_length;
	align.isRC = isRC;
	return align;
}
