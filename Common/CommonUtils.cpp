#include <math.h>
#include "CommonUtils.h"

extDirection oppositeDirection(extDirection dir)
{
	if(dir == SENSE)
	{
		return ANTISENSE;
	}
	else
	{
		return SENSE;
	}	
}

double entropy(const Sequence& s)
{
	std::map<char, int> vals;
	
	int len = s.length();
	for(int i = 0; i < len; i++)
	{
		vals[s[i]]++;
	}
	
	std::map<char, int>::iterator iter;
	
	double sum = 0.0f;
	for(iter = vals.begin(); iter != vals.end(); iter++)
	{
		double f = static_cast<double>(iter->second) / static_cast<double>(len);
		double t = log(f) / log(2.0f);
		
		sum += (f*t);
	}
	
	return -sum;
}

// Make all the permutations of this sequence
void makeExtensions(const PackedSeq& seq, extDirection dir, PSequenceVector& outVector)
{
	
	for(int i = 0; i < NUM_BASES; i++)
	{
		PackedSeq testSeq = seq;
		if(dir == SENSE)
		{
			testSeq.shiftAppend(BASES[i]);
		}
		else
		{
			testSeq.shiftPrepend(BASES[i]);
		}	
		outVector.push_back(testSeq);
	}
}

// Make all the permutations of this sequence
void makePermutations(const Sequence& seq, SequenceVector& outVector)
{
	int seqLen = seq.length();

	// Change every base to a different base
	for(int i = 0; i < seqLen; i++)
	{
		const char calledBase = seq[i];
		for(int j = 0; j < NUM_BASES; j++)
		{
			// copy the sequence to make the permutation
			Sequence perm = seq;
			char newBase = BASES[j];
			if(newBase != calledBase)
			{
				perm[i] = BASES[j];
				outVector.push_back(perm);
			}
		}
	}
} 
