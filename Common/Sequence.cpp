#include "Sequence.h"
#include "Options.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cassert>

using namespace std;

/** Return the complement of the specified base. */
char complementBaseChar(char base)
{
	switch (base) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		default:
			assert(false);
			return 0;
	}
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
	return codeToBase(cstont[baseToCode(anchor)][baseToCode(cs)]);
}

Sequence colourToNucleotideSpace(char anchor, const Sequence& seq)
{
	int seed = baseToCode(anchor);

	stringstream s;
	for (string::const_iterator it = seq.begin();
			it != seq.end(); ++it) {
		seed = cstont[seed][baseToCode(*it)];
		s << codeToBase(seed);
	}
	return s.str();
}
