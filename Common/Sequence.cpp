#include "Sequence.h"

// generate the reverse complement of the sequence
Sequence reverseComplement(const Sequence& s)
{
	Sequence rc;
	
	for(std::string::const_reverse_iterator iter = s.rbegin(); iter != s.rend(); iter++)
	{	
		rc.push_back(complementBaseChar(*iter));
	}
	
	return rc;
}

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
