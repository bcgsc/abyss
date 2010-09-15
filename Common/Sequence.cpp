#include "Sequence.h"
#include "Common/Options.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;

enum { A, C, G, T };
static const int cstont[4][4] = {
	{ A, C, G, T },
	{ C, A, T, G },
	{ G, T, A, C },
	{ T, G, C, A }
};

/** Return the complement of the specified base. */
char complementBaseChar(char base)
{
	switch (base) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		case 'N': return 'N';
		case 'a': return 't';
		case 'c': return 'g';
		case 'g': return 'c';
		case 't': return 'a';
		case 'n': return 'n';
		case '.': return '.';
	}
	cerr << "error: unexpected character: `" << base << "'\n";
	exit(EXIT_FAILURE);
}

/** Return the reverse complement of the specified sequence. */
Sequence reverseComplement(const Sequence& s)
{
	Sequence rc(s);
	reverse(rc.begin(), rc.end());
	if (!opt::colourSpace)
		transform(rc.begin(), rc.end(), rc.begin(),
				complementBaseChar);
	return rc;
}

/** Return the base enumeration for the specified character. */
uint8_t baseToCode(char base)
{
	switch (base) {
		case 'A': case '0': return 0;
		case 'C': case '1': return 1;
		case 'G': case '2': return 2;
		case 'T': case '3': return 3;
	}
	cerr << "error: unexpected character: `" << base << "'\n";
	exit(EXIT_FAILURE);
}

char codeToBase(uint8_t code)
{
	assert(code < 4);
	return (opt::colourSpace ? "0123" : "ACGT")[code];
}

char colourToNucleotideSpace(char anchor, char cs)
{
	return cs == '.' ? 'N'
		: "ACGT"[cstont[baseToCode(anchor)][baseToCode(cs)]];
}

Sequence colourToNucleotideSpace(char anchor, const Sequence& seq)
{
	int seed = baseToCode(anchor);

	ostringstream s;
	s << anchor;
	for (string::const_iterator it = seq.begin();
			it != seq.end(); ++it) {
		seed = cstont[seed][baseToCode(*it)];
		s << codeToBase(seed);
	}
	return s.str();
}

char nucleotideToColourSpace(char a, char b)
{
	if (toupper(a) == 'N' || toupper(b) == 'N')
		return islower(a) || islower(b) ? 'n' : 'N';
	return "0123"[cstont[baseToCode(a)][baseToCode(b)]];
}
