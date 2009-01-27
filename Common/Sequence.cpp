#include <stdio.h>
#include <map>
#include <math.h>
#include "Sequence.h"


// append another string
void seqAppend(Sequence& s, const std::string& str)
{
	s.append(str);
}

// append a base to the string
void seqAppendBase(Sequence& s, const char b)
{
	s.append(1,b);
}

// prepend a base to the string
void seqPrependBase(Sequence& s, const char b)
{
	// stl string doesnt have a prepend function
	std::string temp = s;
	s = b;
	s.append(temp);
	
}

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

// print the sequence
void print(const Sequence& s)
{
	printf("%s", s.c_str());	
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
