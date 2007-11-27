#include <stdio.h>
#include <map>
#include <math.h>
#include "Sequence.h"


// append another string
void seqAppend(Sequence& s, const std::string& str)
{
	s.append(str);
}

// append a base to the string
void seqAppendBase(Sequence& s, const char b)
{
	s.append(1,b);
}

// prepend a base to the string
void seqPrependBase(Sequence& s, const char b)
{
	// stl string doesnt have a prepend function
	std::string temp = s;
	s = b;
	s.append(temp);
	
}

// generate the reverse complement of the sequence
Sequence reverseComplement(const Sequence& s)
{
	Sequence rc;
	
	for(std::string::const_reverse_iterator iter = s.rbegin(); iter != s.rend(); iter++)
	{	
		rc.push_back(complementBase(*iter));
	}
	
	return rc;
}

// print the sequence
void print(const Sequence& s)
{
	printf("%s", s.c_str());	
}


// complement a single base
char complementBase(char base)
{
	switch(base)	
	{
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
	}
	
	printf("could not complement base %c\n", base);
	assert(false && "unknown base");
	return 'N';
}

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
void makeExtensions(const Sequence& seq, extDirection dir, SequenceVector& outVector)
{
	
	int commonLen = seq.length() - 1;
	
	Sequence commonSeq;
	
	if(dir == SENSE)
	{
		commonSeq = seq.substr(1,commonLen);
	}
	else
	{
		commonSeq = seq.substr(0,commonLen);
	}
	
	for(int i = 0; i < NUM_BASES; i++)
	{
		Sequence testSeq = commonSeq;
		if(dir == SENSE)
		{
			seqAppendBase(testSeq, BASES[i]);
		}
		else
		{
			seqPrependBase(testSeq, BASES[i]);
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

