#include "Sequence.h"
#include "Common/Options.h"
#include <algorithm>
#include <cassert>
#include <cctype>
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

/** Return the complement of the specified nucleotide. */
char complementBaseChar(char c)
{
	char rc;
	switch (toupper(c)) {
	  case 'A': rc = 'T'; break;
	  case 'C': rc = 'G'; break;
	  case 'G': rc = 'C'; break;
	  case 'T': rc = 'A'; break;
	  case 'N': rc = 'N'; break;
	  case '.': rc = '.'; break;
	  case 'M': rc = 'K'; break; // A or C
	  case 'R': rc = 'Y'; break; // A or G
	  case 'W': rc = 'W'; break; // A or T
	  case 'S': rc = 'S'; break; // C or G
	  case 'Y': rc = 'R'; break; // C or T
	  case 'K': rc = 'M'; break; // G or T
	  case 'V': rc = 'B'; break; // A or C or G
	  case 'H': rc = 'D'; break; // A or C or T
	  case 'D': rc = 'H'; break; // A or G or T
	  case 'B': rc = 'V'; break; // C or G or T
	  default:
		cerr << "error: unexpected character: `" << c << "'\n";
		exit(EXIT_FAILURE);
	}
	return islower(c) ? tolower(rc) : rc;
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
