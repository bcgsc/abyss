#ifndef COMMONUTILS_H
#define COMMONUTILS_H

#include "Sequence.h"

// A myriad of functions that don't fit well into any class

// make a one base extension of this sequence
void makeExtensions(const Sequence& seq, extDirection dir, SequenceVector& outVector);

// make all the one base permutations of this sequence
void makePermutations(const Sequence& seq, SequenceVector& outVector);

extDirection oppositeDirection(extDirection dir);

// calculate the information entropy of the string
double entropy(const Sequence& s);

#endif
