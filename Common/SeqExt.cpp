#include "SeqExt.h"
#include <cassert>
#include <cstdio>

SeqExt::SeqExt() : m_record(0)
{
}

static uint8_t baseCodeToBit(uint8_t base)
{
	return 1 << base;
}

void SeqExt::setBase(uint8_t base)
{
	m_record |= baseCodeToBit(base);
}

void SeqExt::clearBase(uint8_t base)
{
	m_record &= ~baseCodeToBit(base);
}

bool SeqExt::checkBase(uint8_t base) const
{
	return m_record & baseCodeToBit(base);
}

void SeqExt::ClearAll()
{
	m_record = 0;
}

bool SeqExt::HasExtension() const
{
	return (m_record > 0);
}

bool SeqExt::IsAmbiguous() const
{
	// if the value isn't a power of 2, there is more than one extension
	bool nonZero = m_record > 0;
	bool powerOfTwo = (m_record & (m_record - 1)) > 0;
	return nonZero && powerOfTwo;
}

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
