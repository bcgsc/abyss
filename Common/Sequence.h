#ifndef SEQUENCE_H
#define SEQUENCE_H 1

#include <cstring>
#include <stdint.h>
#include <string>
#include <cassert>

typedef std::string Sequence;

char complementBaseChar(char c);
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
 * Return true if a sequence consists entirely of ACGT chars
 * (case insensitive).
 */
static inline bool allACGT(const Sequence& seq)
{
	return strspn(seq.c_str(), "acgtACGT") == seq.length();
}

/**
 * Transform a sequence in its canonical orientation.
 */
static inline void canonicalize(Sequence& seq)
{
	Sequence rc = reverseComplement(seq);
	if (rc < seq)
		seq = rc;
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

/**
 * Return true if the given sequence contains ambiguity codes.
 */
static inline bool
containsAmbiguityCodes(const Sequence& seq, bool allowN=true)
{
	if (allowN) {
		return seq.find_first_of("MRWSYKVHDBmrwsykvhdbNn")
			!= std::string::npos;
	} else {
		return seq.find_first_of("MRWSYKVHDBmrwsykvhdb")
			!= std::string::npos;
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

/**
 * Overlay one sequence on top of another to create a new sequence.
 * In the cases of differences, the bases in the overlaid sequence
 * take precedence.
 *
 * @param overlay the sequence to be overlaid on target
 * @param target the sequence to be modified/extended
 * @param shift position of overlay sequence relative to target
 * @param maskNew output bases that have been changed or added
 * to target in lowercase.
 */
static inline void overlaySeq(const Sequence& overlay, Sequence& target,
	int shift, bool maskNew = false)
{
	Sequence::const_iterator src = overlay.begin();
	Sequence::iterator dest;

	if (shift < 0) {
		target.insert(0, -shift, 'N');
		dest = target.begin();
	} else {
		assert(shift >= 0);
		if (shift + overlay.length() > target.length()) {
			unsigned suffixLen = shift + overlay.length() -
				target.length();
			target.insert(target.length(), suffixLen, 'N');
		}
		dest = target.begin() + shift;
	}
	for (; src != overlay.end(); ++src, ++dest) {
		assert(dest != target.end());
		if (maskNew && *src != *dest)
			*dest = tolower(*src);
		else
			*dest = *src;
	}
}

#endif
