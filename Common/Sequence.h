#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include <ext/hash_map>
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

// Hash function for Sequences
namespace __gnu_cxx {
template <>
struct hash<Sequence> {
        size_t operator() (const Sequence& x) const {
                return hash<const char*>()(x.c_str());
	// hash<const char*> already exists
        }
};
}

const int SEQUENCE_ID_LENGTH = 32;

char complementBase(char base);
Sequence reverseComplement(const Sequence& s);

// append a base to the string
void seqAppend(Sequence& s, const std::string& str);
void seqAppendBase(Sequence& s, const char b);
void seqPrependBase(Sequence& s, const char b);

#endif
