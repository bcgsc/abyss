#include "ProofReader.h"
#include "PhaseSpace.h"
#include "CommonUtils.h"

#include <fstream>

int main(int argv, char** argc)
{
	Reader fileReader;
	
	if(argv < 2)
	{
		printf("usage: ProofReader <fasta file>\n");
		exit(1);
	}
	
	const char* fastaFile = argc[1];
#if 0
	// Read all the sequences into a vector
	SequenceVector seqVector;
	bool result = fileReader.readFasta(fastaFile, seqVector);
		
	printf("done reading\n");
		
	// assume reads are all the same length
	int readLen = seqVector.front().length();

	// record all the multiplicity of sequences
	SeqRecord multiplicity;
	multiplicity.addMultipleSequences(seqVector);
	
	printf("done multiplicity load\n");
	
	// map of sequence->correction
	std::map<Sequence, Sequence> corrections;

	correctReads(multiplicity, corrections);
	printf("done correction\n");
	
	outputCorrectedSequences(seqVector, multiplicity, corrections);
	printf("done output\n");
#endif
	return 0;
}

// Correct the reads by finding sequences with multiplicity 1. If the read has a multiplicity 1 and a unique, non-singular correction, make the correction
// The correction is performed in place
void correctReads(const SeqRecord& seqMult, std::map<Sequence, Sequence>& corrections)
{
	// perform the actual correction
	int numElems = seqMult.size();
	int count = 0;
	
	for(ConstSeqRecordIter iter = seqMult.begin(); iter != seqMult.end(); iter++)
	{
		if(count % 10000 == 0)
		{
			printf("correction: %d/%d\n", count, numElems);
		}
		
		// is the read multiplicity 1?
		if(iter->second == 1)
		{
			// this sequence is potentially correctable
			Sequence corrResult = correctSequence(iter->first, seqMult);
			if(corrResult != iter->first)
			{
				// sequence was corrected
				corrections[iter->first] = corrResult;
			}
		}

		count++;
	}
}

// take in a sequence and a prb file and correct it if possible, otherwise return the original sequence
Sequence correctSequence(const Sequence& seq, const SeqRecord& multiplicity)
{
	int seqLen = seq.length();
	
	int numCorrections = 0;
	Sequence correction;
	
	SequenceVector perms;
	makePermutations(seq, perms);
	
	for(SequenceVector::const_iterator iter = perms.begin(); iter != perms.end(); iter++)
	{
		if(multiplicity.getMultiplicity(*iter) > 1)
		{
			numCorrections++;
			correction = *iter;
		}
	}
	
	if(numCorrections == 1)
	{
		return correction;
	}
	else
	{
		return seq;
	}
}

// this function outputs the reads for assembly
void outputCorrectedSequences(const SequenceVector& seqVector, const SeqRecord& multiplicity, std::map<Sequence, Sequence>& corrections)
{			
	std::fstream outFasta;
	outFasta.open("autocorrect.fa", std::ios::out);
	
	int outputCount = 0;
	int count = 0;
	int numReads = seqVector.size();
	for(ConstSequenceVectorIterator iter = seqVector.begin(); iter != seqVector.end(); iter++)
	{
		Sequence s = *iter;
		if(corrections.find(s) != corrections.end())
		{
			//sequence has correction
			s = corrections.find(s)->second;
		}
		
		int mult = multiplicity.getMultiplicity(s);
		
		outFasta << ">" << outputCount << " " << mult << "\n" << s << "\n";
		outputCount++;
	}
	outFasta.close();
}

// Look for high-quality sequences that prematurely end with no extension. Attempt to fill that hole with a low-quality sequence that can be permutated
void fillHoles(const SequenceVector& seqVector, const SeqRecord& multiplicity, const PhaseSpace& phase, std::map<Sequence, Sequence>& corrections)
{
	assert(seqVector.size() > 0);
	
	// perform the actual correction
	int numElems = seqVector.size();
	
	int count = 0;
	for(int i = 0; i < numElems; i++)
	{
		if(i % 1000 == 0)
		{
			//printf("%d\n", i);
		}
		
		// get the sequence
		const Sequence& seq = seqVector[i];
		//const ReadPrb& readPrb = prbVector[i];

		// Check that this sequence is reasonably high quality
		//double Perror = calculatePError(readPrb);
		
		//if(Perror < 0.1f)
		{
			// this sequence is high quality
			bool hasParent = phase.hasParent(seq);
			bool hasChild = phase.hasChild(seq);
			
			if((hasParent && hasChild) || (!hasParent && !hasChild))
			{
				continue;
			}
			else
			{
				extDirection dir;	
				if(!hasParent)
				{
					dir = ANTISENSE;
				}
				else
				{
					dir = SENSE;
				}
				
				// create the possible extensions
				SequenceVector extensions;
				makeExtensions(seq, dir, extensions);
				
				for(ConstSequenceVectorIterator iter = extensions.begin(); iter != extensions.end(); iter++)
				{
					SequenceVector permutations;
					makePermutations(*iter, permutations);
					
					int numCorrections = 0;
					for(SequenceVector::const_iterator iter2 = permutations.begin(); iter2 != permutations.end(); iter2++)
					{
						// check if this is a good correction (the permutation exists, has mult = 1 and no parents/children
						if(multiplicity.getMultiplicity(*iter2) == 1)
						{
							//printf("correcting read\n");
							corrections[*iter2] = *iter;
						}
					}
				}
			}	
		}	
		count++;
	}
}

