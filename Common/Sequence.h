#ifndef SEQUENCE_H
#define SEQUENCE_H 1

#include <string>
#include <vector>

typedef std::string Sequence;

Sequence reverseComplement(const Sequence& s);
Sequence colourToNucleotideSpace(char anchor, const Sequence& seq);

// Create the two bit code for the base
uint8_t baseToCode(char base);
char codeToBase(uint8_t code);

typedef std::vector<Sequence> SequenceVector;
typedef SequenceVector::iterator SequenceVectorIterator;

enum { A, C, G, T };
static const int cstont[4][4] = {
	{ A, C, G, T },
	{ C, A, T, G },
	{ G, T, A, C },
	{ T, G, C, A }
};

#endif
