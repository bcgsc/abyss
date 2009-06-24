#include "Sequence.h"
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
	transform(rc.begin(), rc.end(), rc.begin(), complementBaseChar);
	return rc;
}
