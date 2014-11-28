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

/** Return true if c is one of [ACGTacgt]. */
static inline bool isACGT(char c)
{
	return c == 'A' || c == 'C' || c == 'G' || c == 'T'
		|| c == 'a' || c == 'c' || c == 'g' || c == 't';
}

/**
 * Convert each ambiguity code to the lexicographically smallest
 * matching base.
 */
static inline void
flattenAmbiguityCodes(Sequence& seq, bool skipNs=true)
{
	for (Sequence::iterator it = seq.begin(); it != seq.end(); ++it) {
		switch (toupper(*it)) {
			case 'N':
				if (!skipNs)
					*it = 'A';
				break;
			case 'M': *it = 'A'; break;
			case 'R': *it = 'A'; break;
			case 'W': *it = 'A'; break;
			case 'S': *it = 'C'; break;
			case 'Y': *it = 'C'; break;
			case 'K': *it = 'G'; break;
			case 'V': *it = 'A'; break;
			case 'H': *it = 'A'; break;
			case 'D': *it = 'A'; break;
			case 'B': *it = 'C'; break;
			default:
				break;
		}
	}
}

unsigned ambiguityToBitmask(char c);
unsigned bitmaskToAmbiguity(unsigned x);

/** Return the bitwise-and of the specified ambiguity codes. */
static inline char ambiguityAnd(char ca, char cb)
{
	char c = bitmaskToAmbiguity(
			ambiguityToBitmask(ca) & ambiguityToBitmask(cb));
	return islower(ca) && islower(cb) ? tolower(c) : c;
}

/** Return the bitwise-or of the specified ambiguity codes. */
static inline char ambiguityOr(char ca, char cb)
{
	char c = bitmaskToAmbiguity(
			ambiguityToBitmask(ca) | ambiguityToBitmask(cb));
	return islower(ca) || islower(cb) ? tolower(c) : c;
}

/** Return whether one ambiguity code is a subset of the other.
 */
static inline bool ambiguityIsSubset(char a, char b)
{
	char intersection = ambiguityAnd(a, b);
	return intersection == a || intersection == b;
}

#endif
