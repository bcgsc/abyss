#include "SeqExt.h"
#include <cassert>
#include <cstdio>

SeqExt SeqExt::complement() const
{
	SeqExt comp;
	for (int i = 0; i < NUM_BASES; i++)
		if (checkBase(i))
			comp.setBase(complementBaseCode(i));
	return comp;
}

void SeqExt::print() const
{
	assert(m_record < 1<<NUM_BASES);
	printf("ext: %c%c%c%c\n",
			 checkBase(3) ? 'T' : ' ',
			 checkBase(2) ? 'G' : ' ',
			 checkBase(1) ? 'C' : ' ',
			 checkBase(0) ? 'A' : ' ');
}
