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

/** Construct a k-mer pair from two k-mer. */
KmerPair(const Kmer& a, const Kmer& b) : m_a(a), m_b(b) { }

/** Construct a k-mer pair from two strings. */
KmerPair(const std::string& a, const std::string& b) : m_a(a), m_b(b) { }

/** Construct a k-mer pair from one string.
 * The first and last word of the specified string are used to construct the
 * two k-mers. The two words may overlap.
 */
KmerPair(const std::string& s) :
	m_a(s.substr(0, Kmer::length())),
	m_b(s.substr(s.size() - Kmer::length(), Kmer::length()))
{
}

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
static unsigned length()
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
	return m_a == ::reverseComplement(m_b);
}

/** Return whether the specified k-mer pair edge is palindromic. */
bool isPalindrome(extDirection dir) const
{
	Kmer a(m_a);
	if (dir == SENSE)
		a.shift(SENSE, 0);
	else
		a.setLastBase(SENSE, 0);

	Kmer b(m_b);
	b.reverseComplement();
	if (dir == ANTISENSE)
		b.shift(SENSE, 0);
	else
		b.setLastBase(SENSE, 0);

	return a == b;
}

/** Return a string representation of the k-mer pair, using an accurate number
 * of Ns to separate the two k-mer.
 */
std::string str() const
{
	assert(length() >= m_a.length());
	std::string s(length(), 'N');
	s.replace(0, m_a.length(), m_a.str());
	s.replace(length() - m_b.length(), m_b.length(), m_b.str());
	return s;
}

/** Return a string representation of the k-mer pair, using the specified
 * separator string to separate the two k-mer.
 */
std::string str(const char* sep) const
{
	std::string s;
	s.reserve(2 * Kmer::length() + 1);
	s += m_a.str();
	s += sep;
	s += m_b.str();
	return s;
}

/** Set the last base of each k-mer. */
void setLastBase(extDirection sense, const Dinuc& x)
{
	m_a.setLastBase(sense, x.a());
	m_b.setLastBase(sense, x.b());
}

/** Shift both k-mer. */
Dinuc shift(extDirection sense, Dinuc x = Dinuc(0))
{
	return Dinuc(
			m_a.shift(sense, x.a()),
			m_b.shift(sense, x.b()));
}

/** Reverse complement this k-mer pair. */
void reverseComplement()
{
	m_a.reverseComplement();
	m_b.reverseComplement();
	std::swap(m_a, m_b);
}

/** Print this k-mer pair. */
friend std::ostream& operator<<(std::ostream& out, const KmerPair& x)
{
	return out << x.m_a.str() << '-' << x.m_b.str();
}

/** Return the number of bytes needed. */
static unsigned serialSize()
{
	return sizeof(KmerPair);
}

/** Serialize this k-mer pair. */
size_t serialize(void* dest) const
{
	memcpy(dest, this, sizeof *this);
	return sizeof *this;
}

/** Unserialize this k-mer pair. */
size_t unserialize(const void* src)
{
	memcpy(this, src, sizeof *this);
	return sizeof *this;
}

/** Return a hash value that does not change with reverse complementation. */
unsigned getCode() const
{
	return m_a.getCode() ^ m_b.getCode();
}

private:

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
