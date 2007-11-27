#include "ProofReader.h"
#include "PhaseSpace.h"


int main(int argv, char** argc)
{
	Reader fileReader;
	
	if(argv < 2)
	{
		printf("usage: ProofReader <fasta file> <apb file>\n");
		exit(1);
	}
	
	const char* fastaFile = argc[1];
	
	// apb is an annotated prb file, simply a PRB that has been processed with seperateReads.pl to split into two and give the same header as the fasta seq file
	// the seq and apb file must be a line for line match or else really bad things will happen TODO: check for this
	const char* apbFile = argc[2];
	
	// Read all the sequences into a vector
	SequenceVector seqVector;
	bool result = fileReader.readFasta(fastaFile, seqVector);
	
	// Read all the prb records into a vector
	PrbVector prbVector;
	prbVector.reserve(seqVector.size());
	//result = fileReader.readAPB(apbFile, prbVector);
	
	// record all the multiplicity of sequences
	SeqRecord multiplicity;
	multiplicity.addMultipleSequences(seqVector);
	
	// map of original sequence -> corrected sequence
	std::map<Sequence, Sequence> corrections;	

	correctReads(seqVector, prbVector, multiplicity);
	//fillHoles(seqVector, prbVector, multiplicity, corrections);
	outputReadSet(seqVector, prbVector, multiplicity, corrections);
	return 0;
}

// Correct the reads by finding low-quality sequences with multiplicity 1. If the read has a multiplicity 1 and a unique, non-singular correction, make the correction
// The correction is performed in place
void correctReads(SequenceVector& seqVector, const PrbVector& prbVector, const SeqRecord& multiplicity)
{
	// perform the actual correction
	int numElems = seqVector.size();
	int count = 0;
	for(int i = 0; i < numElems; i++)
	{
		// get the sequence
		const Sequence& seq = seqVector[i];
		const ReadPrb& readPrb = prbVector[i];

		if(multiplicity.getMultiplicity(seq) == 1)
		{
			// this sequence is potentially correctable
			seqVector[i] = correctSequence(seq, readPrb, multiplicity);			
		}
	}
}

// this function outputs the reads for assembly
void outputReadSet(const SequenceVector& seqVector, const PrbVector& prbVector, const SeqRecord& multiplicity, std::map<Sequence, Sequence>& corrections)
{
	std::map<Sequence, bool> seenSeq;
	
	// assume reads are all the same length
	int readLen = seqVector.front().length();
		
	// should cache this
	PhaseSpace phase(readLen);
	phase.addReads(seqVector);
			
	int count = 0;
	for(ConstSequenceVectorIterator iter = seqVector.begin(); iter != seqVector.end(); iter++)
	{
		Sequence s = *iter;
		if(corrections.count(*iter) > 0)
		{
			s = corrections[*iter];
		}
		
		bool hasParent = phase.hasParent(s);
		bool hasChild = phase.hasChild(s);		
		if(hasParent && hasChild && !seenSeq.count(s) > 0)
		{
			printf(">%d\n%s\n", count, s.c_str());
			count++;
			
			seenSeq[s] = 1;
		}
	}
}

// Look for high-quality sequences that prematurely end with no extension. Attempt to fill that hole with a low-quality sequence that can be permutated
void fillHoles(const SequenceVector& seqVector, const PrbVector& prbVector, const SeqRecord& multiplicity, std::map<Sequence, Sequence>& corrections)
{
	assert(seqVector.size() > 0);
	
	// assume reads are all the same length
	int readLen = seqVector.front().length();
	
	// Create the phase space
	PhaseSpace phase(readLen);
	
	phase.addReads(seqVector);
	
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

// take in a sequence and a prb file and correct it if possible, otherwise return the original sequence
Sequence correctSequence(const Sequence& seq, const ReadPrb& readPrb, const SeqRecord& multiplicity)
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

