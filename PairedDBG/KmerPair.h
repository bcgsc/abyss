#ifndef PAIREDDBG_KMERPAIR_H
#define PAIREDDBG_KMERPAIR_H 1

#include <stdint.h>

#include "Dinuc.h"
#include "Common/Kmer.h"

/** A pair of k-mer. */
class KmerPair
{
public:
	typedef Dinuc::Nuc Nuc;

/** Return whether the two objects are equal. */
bool operator==(const KmerPair& x) const
{
	return m_a == x.m_a && m_b == x.m_b;
}

/** Return whether the two objects are inequal. */
bool operator!=(const KmerPair& x) const
{
	return !(*this == x);
}

/** Return whether this object is less than the other. */
bool operator<(const KmerPair& x) const
{
	return m_a != x.m_a
		? m_a < x.m_a
		: m_b < x.m_b;
}

/** Return the length of this k-mer pair. */
unsigned length() const
{
	//xxx fixme This is not correct.
	return m_a.length() + m_b.length();
}

/** Return the terminal nucleotides. */
char getLastBaseChar() const
{
	//xxx fixme This is not correct.
	return m_a.getLastBaseChar();
}

/** Return the hash value. */
uint64_t getHashCode() const
{
	return m_a.getHashCode() ^ m_b.getHashCode();
}

/** Return whether this k-mer pair is palindromic. */
bool isPalindrome() const
{
	//xxx fixme Is this correct?
	return m_a.isPalindrome() && m_b.isPalindrome();
}

/** Return a string represenation. */
std::string str() const
{
	//xxx fixme todo
	return std::string("xxx");
}

/** Set the last base of each k-mer. */
void setLastBase(extDirection sense, const Dinuc& x)
{
	m_a.setLastBase(sense, x.a());
	m_b.setLastBase(sense, x.b());
}

/** Shift both k-mer. */
Dinuc shift(extDirection sense, Dinuc x = 0)
{
	return sense == SENSE ? shiftAppend(x) : shiftPrepend(x);
}

/** Reverse complement this k-mer pair. */
void reverseComplement()
{
	m_a.reverseComplement();
	m_b.reverseComplement();
	std::swap(m_a, m_b); //xxx Is this correct?
}

private:

/** Shift both k-mer left. */
Dinuc shiftAppend(Dinuc x)
{
	Nuc a = m_a.shift(SENSE, x.a());
	Nuc b = m_b.shift(SENSE, x.b());
	return Dinuc(a, b);
}

/** Shift both k-mer right. */
Dinuc shiftPrepend(Dinuc x)
{
	Nuc a = m_a.shift(ANTISENSE, x.a());
	Nuc b = m_b.shift(ANTISENSE, x.b());
	return Dinuc(a, b);
}

	/** The first k-mer. */
	Kmer m_a;

	/** The second k-mer. */
	Kmer m_b;
};


/** Return the reverse complement. */
static inline KmerPair reverseComplement(const KmerPair& u)
{
	KmerPair urc(u);
	urc.reverseComplement();
	return urc;
}

NAMESPACE_STD_HASH_BEGIN
	template <> struct hash<KmerPair> {
		size_t operator()(const KmerPair& u) const
		{
			return u.getHashCode();
		}
	};
NAMESPACE_STD_HASH_END

#endif
