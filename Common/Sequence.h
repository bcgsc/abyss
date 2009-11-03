#ifndef SEQUENCE_H
#define SEQUENCE_H 1

#include <stdint.h>
#include <string>

typedef std::string Sequence;

Sequence reverseComplement(const Sequence& s);
Sequence colourToNucleotideSpace(char anchor, const Sequence& seq);
char colourToNucleotideSpace(char anchor, char cs);
char nucleotideToColourSpace(char a, char b);

uint8_t baseToCode(char base);
char codeToBase(uint8_t code);

#endif
