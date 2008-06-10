#include "SeqExt.h"

SeqExt::SeqExt() : m_record(0)
{
	
}

//
//
//
void SeqExt::SetBase(char base)
{
	unsigned char bit = base2Bit(base);
	m_record |= bit;
}

//
//
//
void SeqExt::ClearBase(char base)
{
	unsigned char bit = base2Bit(base);
	unsigned char mask = ~bit;
	m_record &= mask;
}

//
//
//
bool SeqExt::CheckBase(char base) const
{
	unsigned char bit = base2Bit(base);
	return (m_record & bit);
}

//
//
//
void SeqExt::ClearAll()
{
	m_record = 0;
}

//
//
//
bool SeqExt::HasExtension() const
{
	return (m_record > 0);
}

//
//
//
bool SeqExt::IsAmbiguous() const
{
	// if the value isn't a power of 2, there is more than one extension
	bool nonZero = m_record > 0;
	bool powerOfTwo = (m_record & (m_record - 1)) > 0;
	return nonZero && powerOfTwo;
}

//
//
//
SeqExt SeqExt::complement() const
{
	SeqExt comp;
	for(int i = 0; i < NUM_BASES; i++)
	{
		char currBase = BASES[i];
		if(CheckBase(currBase))
		{
			comp.SetBase(::complementBaseChar(currBase));
		}
	}
	return comp;
}

//
void SeqExt::print() const
{
	printf("ext: %c%c%c%c\n",
			 CheckBase('T') ? 'T' : ' ',
			 CheckBase('G') ? 'G' : ' ',
			 CheckBase('C') ? 'C' : ' ',
			 CheckBase('A') ? 'A' : ' ');
}

//
//
//
unsigned char SeqExt::base2Bit(char base)
{
	if(base == 'A')
	{
		return 0x1;
	}
	else if(base == 'C')
	{
		return 0x2;
	}
	else if(base == 'G')
	{
		return 0x4;
	}
	else if(base == 'T')
	{
		return 0x8;
	}
	else
	{
		assert(false);	
		return 0x1;
	}
}

//
//
//
char SeqExt::bit2Base(unsigned char code)
{
	if(code == 0x1)
	{
		return 'A';
	}
	else if(code == 0x2)
	{
		return 'C';
	}
	else if(code == 0x4)
	{
		return 'G';
	}
	else if(code == 0x8)
	{
		return 'T';
	}
	else
	{
		// unknown code
		assert(false);
		return 'A';
	}
}
