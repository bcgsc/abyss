#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include "CommonDefs.h"

enum SequenceAdjacency
{
	SA_NONE,
	SA_SAME_SENSE,
	SA_SAME_ANTISENSE,
	SA_RC_SENSE,
	SA_RC_ANTISENSE	
};

// SENSE AND ANTISENSE HAVE TO BE ZERO AND ONE
enum extDirection
{
	SENSE = 0,
	ANTISENSE = 1,
	NUM_DIRECTIONS
};

const int SEQUENCE_ID_LENGTH = 32;

char complementBase(char base);
Sequence reverseComplement(const Sequence& s);

// append a base to the string
void seqAppend(Sequence& s, const std::string& str);
void seqAppendBase(Sequence& s, const char b);
void seqPrependBase(Sequence& s, const char b);

// make a one base extension of this sequence
void makeExtensions(const Sequence& seq, extDirection dir, SequenceVector& outVector);

// make all the one base permutations of this sequence
void makePermutations(const Sequence& seq, SequenceVector& outVector);

extDirection oppositeDirection(extDirection dir);

// calculate the information entropy of the string
double entropy(const Sequence& s);

#endif
