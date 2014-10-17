#ifndef PAIREDDBG_DINUC_H
#define PAIREDDBG_DINUC_H 1

#include "Common/BitUtil.h"
#include <cassert>
#include <stdint.h>

/** A pair of nucleotides. */
class Dinuc
{
public:
	/** The number of symbols. */
	static const unsigned NUM = 16;

	/** A nucleotide. A bit vector of two bits. */
	typedef uint8_t Nuc;

	/** A dinucleotide. A bit vector of four bits. */
	typedef uint8_t Bits;

	/** Construct a Dinuc from two nucleotides. */
	Dinuc(Nuc a, Nuc b) : m_data(a | (b << 2)) { }

	/** Construct a Dinuc from an integer. */
	explicit Dinuc(Bits x) : m_data(x) { }

	/** Cast to an integer. */
	Bits toInt() const { return m_data; }

	/** Return the first nucleotide. */
	Nuc a() const { return m_data & 0x3; }

	/** Return the first nucleotide. */
	Nuc b() const { return (m_data >> 2) & 0x3; }

	/** Compare two dinucleotides. */
	bool operator<(const Dinuc& x) const
	{
		return m_data < x.m_data;
	}

	/** Complement a single base. */
	static Nuc complementNuc(Nuc x) { return 3 - x; }

	/** Return the reverse complement of this dinucleotide. */
	Dinuc reverseComplement() const
	{
		return Dinuc(complementNuc(b()), complementNuc(a()));
	}

	/** Increment this dinucleotide. */
	Dinuc& operator++()
	{
		++m_data;
		return *this;
	}

	/** Return the first dinucleotide. */
	static Dinuc begin() { return Dinuc(0); }

	/** Return the last dinucleotide. */
	static Dinuc end() { return Dinuc(NUM); }

private:
	/** Two nucleotides packed into a single scalar. */
	Bits m_data;
};

/** Return the reverse complement of this dinucleotide. */
static inline Dinuc reverseComplement(const Dinuc& x)
{
	return x.reverseComplement();
}

/** A set of dinucleotides. */
class DinucSet
{
public:
	typedef Dinuc Symbol;

	/** The number of symbols. */
	static const unsigned NUM = Dinuc::NUM;

	/** A bit vector. */
	typedef uint16_t Bits;

	/** Default constructor. */
	DinucSet() : m_data(0) { }

	/** Construct a set containing a single element. */
	DinucSet(const Dinuc& x) : m_data(1 << x.toInt()) { }

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
	return m_data & (1 << x.toInt());
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
	m_data |= 1 << x.toInt();
}

/** Remove all elements from this set. */
void clear()
{
	m_data = 0;
}

/** Remove the specified elements from this set. */
void clear(const DinucSet& x)
{
	m_data &= ~x.m_data;
}

/** Return the complementary nucleotides of this set. */
DinucSet complement() const
{
	DinucSet x;
	for (Dinuc i = Dinuc::begin(); i < Dinuc::end(); ++i) {
		if (checkBase(i))
			x.setBase(i.reverseComplement());
	}
	return x;
}

bool operator==(const DinucSet& x) const
{
	return m_data == x.m_data;
}

private:
	/** A bit vector representing a set. */
	Bits m_data;
};

#endif
