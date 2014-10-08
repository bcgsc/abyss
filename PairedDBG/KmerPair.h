#ifndef PAIREDDBG_KMERPAIR_H
#define PAIREDDBG_KMERPAIR_H 1

#include "Dinuc.h"

#include "Common/Kmer.h"

#include <stdint.h>
#include <utility>

/** A pair of k-mer. */
class KmerPair
{
public:
	typedef Dinuc::Nuc Nuc;

/** Default constructor. */
KmerPair() { }

/** Construct a k-mer pair. */
KmerPair(const Kmer& a, const Kmer& b) : m_a(a), m_b(b) { }

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

/** Return the length of a the k-mer pair, including the gap. */
unsigned length() const
{
	return s_length;
}

/** Set the length of a k-mer pair, including the gap.
 * This value is shared by all instances.
 */
static void setLength(unsigned length)
{
	assert(length >= 2 * Kmer::length());
	s_length = length;
}

/** Return the terminal nucleotides. */
std::pair<char, char> getLastBaseChar() const
{
	return std::make_pair(m_a.getLastBaseChar(), m_b.getLastBaseChar());
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

/** Return whether the specified k-mer pair edge is palindromic. */
bool isPalindrome(extDirection dir) const
{
	//xxx fixme Is this correct?
	return m_a.isPalindrome(dir) && m_b.isPalindrome(dir);
}

/** Return a string represenation. */
std::string str() const
{
	assert(length() >= m_a.length());
	std::string s(length(), 'N');
	s.replace(0, m_a.length(), m_a.str());
	s.replace(length() - m_b.length(), m_b.length(), m_b.str());
	return s;
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

friend std::ostream& operator<<(std::ostream& out, const KmerPair& x)
{
	return out << x.m_a.str() << '-' << x.m_b.str();
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

	/** The length of a k-mer pair, including the gap. */
	static unsigned s_length;

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
