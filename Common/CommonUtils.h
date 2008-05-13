#ifndef COMMONUTILS_H
#define COMMONUTILS_H

#include "Sequence.h"

// A myriad of functions that don't fit well into any class
extDirection oppositeDirection(extDirection dir);

// calculate the information entropy of the string
double entropy(const Sequence& s);

// min/max
int min(const int& n1, const int& n2);
int max(const int& n1, const int& n2);

// complement a base
char complementBaseChar(char b);

static inline uint8_t complementBaseCode(uint8_t b)
{
	return ~b & 0x3;
}

void PrintBufferAsHex(char* buffer, int length);

#endif
