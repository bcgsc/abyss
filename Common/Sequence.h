#ifndef SEQUENCE_H
#define SEQUENCE_H 1

#include <string>
#include <vector>

typedef std::string Sequence;

Sequence reverseComplement(const Sequence& s);

typedef std::vector<Sequence> SequenceVector;
typedef SequenceVector::iterator SequenceVectorIterator;

#endif
