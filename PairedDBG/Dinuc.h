#ifndef PAIREDDBG_DINUC_H
#define PAIREDDBG_DINUC_H 1

#include "Common/BitUtil.h"
#include <cassert>
#include <stdint.h>

/** A pair of nucleotides. */
class Dinuc
{
public:
	/** A nucleotide. A bit vector of two bits. */
	typedef uint8_t Nuc;

	/** A dinucleotide. A bit vector of four bits. */
	typedef uint8_t Bits;

	/** Construct a Dinuc from two nucleotides. */
	Dinuc(Nuc a, Nuc b) : m_data(a | (b << 2)) { }

	/** Construct a Dinuc from an integer. */
	Dinuc(Bits x) : m_data(x) { }

	/** Cast to an integer. */
	operator Bits() const { return m_data; }

	/** Return the first nucleotide. */
	Nuc a() const { return m_data & 0x3; }

	/** Return the first nucleotide. */
	Nuc b() const { return (m_data >> 2) & 0x3; }

private:
	/** Two nucleotides packed into a single scalar. */
	Bits m_data;
};

/** A set of dinucleotides. */
class DinucSet
{
public:
	/** The maximum number of out edges. */
	static const unsigned NUM_EDGES = 16;

	/** A bit vector. */
	typedef uint16_t Bits;

	DinucSet() : m_data(0) { }
	DinucSet(const Dinuc& x) : m_data(1 << x) { }

/** Return a set with the specified bits set. */
static DinucSet mask(Bits x)
{
	DinucSet s;
	s.m_data = x;
	return s;
}

/** Return whether the specified element is present in this set. */
bool checkBase(const Dinuc& x) const
{
	return m_data & (1 << x);
}

/** Return the number of elements in this set. */
unsigned outDegree() const
{
	return popcount(m_data);
}

/** Return whether this set is non-empty. */
bool hasExtension() const
{
	return m_data != 0;
}

/** Return whether this set has two or more elements. */
bool isAmbiguous() const
{
	return outDegree() > 1;
}

/** Add the specified element to this set. */
void setBase(const Dinuc& x)
{
	m_data |= 1 << x;
}

/** Remove the specified elements from this set. */
void clear(const DinucSet& x)
{
	m_data &= ~x.m_data;
}

/** Return the complementary nucleotides of this set. */
DinucSet complement() const
{
	//xxx fixme todo
//	assert(false);
	uint16_t val=this->Bits;
	DinucSet Ds=mask(reverseBits(val));
	return Ds;
}
//can be improved
uint16_t reverseBits(uint16_t val) const
{
	uint16_t y=0;
	 for(int position=NUM_DINU-1; position>=0; position--){
	y+=((val&1)<<position);
	val >>= 1;
	assert(position < 1<<NUM_DINU);
	}	 
	 return y;
}
private:
	/** A bit vector representing a set. */
	Bits m_data;
};

#endif
