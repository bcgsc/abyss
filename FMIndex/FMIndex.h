#ifndef FMINDEX_H
#define FMINDEX_H 1

#include "config.h"
#include "BitArrays.h"
#include "IOUtil.h"
#include "sais.hxx"
#include <boost/integer.hpp>
#include <algorithm>
#include <cassert>
#include <cstdlib> // for exit
#include <iostream>
#include <iterator>
#include <limits> // for numeric_limits
#include <stdint.h>
#include <string>
#include <utility>
#include <vector>

/** A match of a substring of a query sequence to an FM index. */
struct FMInterval
{
	FMInterval(size_t l, size_t u, unsigned qstart, unsigned qend)
		: l(l), u(u), qstart(qstart), qend(qend) { }
	size_t l, u;
	unsigned qstart, qend;

	unsigned qspan() const
	{
		assert(qstart <= qend);
		return qend - qstart;
	}
};

/** A match of a substring of a query sequence to a target sequence.
 */
struct Match
{
	unsigned qstart, qend;
	size_t tstart;
	unsigned count;
	Match(unsigned qstart, unsigned qend,
			size_t tstart, unsigned count)
		: qstart(qstart), qend(qend), tstart(tstart), count(count) { }

	unsigned qspan() const
	{
		assert(qstart <= qend);
		return qend - qstart;
	}
};

/** An FM index. */
class FMIndex
{
	/** A symbol. */
	typedef uint8_t T;

	/** The sentinel symbol. */
	static T SENTINEL() { return std::numeric_limits<T>::max(); }

  public:
	/** An index. */
	typedef boost::uint_t<FMBITS>::least size_type;

	/** An index for SAIS, which must be signed. */
	typedef boost::int_t<FMBITS>::least sais_size_type;

	/** The type of a symbol. */
	typedef T value_type;

FMIndex() : m_sampleSA(1) { }

/** Return the size of the string not counting the sentinel. */
size_t size() const { return m_occ.size() - 1; }

/** Build an FM-index of the specified file. */
template<typename It>
void assign(It first, It last)
{
	assert(first < last);
	assert(size_t(last - first)
			< std::numeric_limits<size_type>::max());

	m_sampleSA = 1;

	// Translate the alphabet.
	if (m_alphabet.empty())
		setAlphabet(first, last);
	std::transform(first, last, first, Translate(*this));
	std::replace(first, last, SENTINEL(), T(0));

	std::cerr << "The alphabet has "
		<< m_alphabet.size() << " symbols.\n";

	// Construct the suffix array.
	std::cerr << "Building the suffix array...\n";
	size_t n = last - first;
	m_sa.resize(n + 1);
	m_sa[0] = n;

	assert(sizeof (size_type) == sizeof (sais_size_type));
	int status = saisxx(first,
			reinterpret_cast<sais_size_type*>(&m_sa[1]),
			(sais_size_type)n,
			(sais_size_type)m_alphabet.size());
	assert(status == 0);
	if (status != 0)
		abort();

	// Construct the Burrows-Wheeler transform.
	std::vector<T> bwt;
	std::cerr << "Building the Burrows-Wheeler transform...\n";
	bwt.resize(m_sa.size());
	for (size_t i = 0; i < m_sa.size(); i++)
		bwt[i] = m_sa[i] == 0 ? SENTINEL() : first[m_sa[i] - 1];

	std::cerr << "Building the character occurence table...\n";
	m_occ.assign(bwt);
	countOccurences();
}

/** Sample the suffix array. */
void sampleSA(unsigned period)
{
	assert(period > 0);
	if (period == m_sampleSA)
		return;
	assert(m_sampleSA == 1);
	m_sampleSA = period;
	if (m_sampleSA == 1)
		return;
	std::vector<size_type>::iterator out = m_sa.begin();
	for (size_t i = 0; i < m_sa.size(); i += m_sampleSA)
		*out++ = m_sa[i];
	m_sa.erase(out, m_sa.end());
	assert(!m_sa.empty());
}

/** Return the position of the specified suffix in the original
 * string.
 */
size_t locate(size_t i) const
{
	size_t n = 0;
	while (i % m_sampleSA != 0) {
		T c = m_occ.at(i);
		i = c == SENTINEL() ? 0 : m_cf[c] + m_occ.rank(c, i);
		n++;
	}
	assert(i / m_sampleSA < m_sa.size());
	size_t pos = m_sa[i / m_sampleSA] + n;
	return pos < m_occ.size() ? pos : pos - m_occ.size();
}

/** Decompress the index. */
template <typename It>
void decompress(It out)
{
	// Construct the original string.
	std::vector<T> s;
	for (size_t i = 0;;) {
		assert(i < m_occ.size());
		T c = m_occ.at(i);
		if (c == SENTINEL())
			break;
		s.push_back(c);
		i = m_cf[c] + m_occ.rank(c, i);
		assert(i > 0);
	}

	// Translate the character set and output the result.
	for (std::vector<T>::const_reverse_iterator it = s.rbegin();
			it != s.rend(); ++it) {
		T c = *it;
		assert(c < m_alphabet.size());
		*out++ = m_alphabet[c];
	}
}

/** Search for an exact match. */
template <typename It>
std::pair<size_t, size_t> findExact(It first, It last) const
{
	assert(first < last);
	size_t l = 1, u = m_occ.size();
	It it;
	for (it = last - 1; it >= first && l < u; --it) {
		T c = *it;
		if (c == SENTINEL())
			return std::make_pair(0, 0);
		l = m_cf[c] + m_occ.rank(c, l);
		u = m_cf[c] + m_occ.rank(c, u);
	}
	return std::make_pair(l, u);
}

/** Search for an exact match. */
template <typename String>
std::pair<size_t, size_t> findExact(const String& q) const
{
	String s(q.size());
	std::transform(q.begin(), q.end(), s.begin(), Translate(*this));
	return findExact(s.begin(), s.end());
}

/** Search for a suffix of the query that matches a prefix of the
 * target.
 */
template <typename It, typename OutIt>
OutIt findOverlapSuffix(It first, It last, OutIt out,
		unsigned minOverlap) const
{
	assert(first < last);

	size_t l = 1, u = m_occ.size();
	for (It it = last - 1; it >= first && l < u; --it) {
		T c = *it;
		if (c == SENTINEL())
			break;
		l = m_cf[c] + m_occ.rank(c, l);
		u = m_cf[c] + m_occ.rank(c, u);
		if (l >= u)
			break;

		if (last - it < minOverlap)
			continue;

		const char sep = 0;
		size_t l1 = m_cf[sep] + m_occ.rank(sep, l);
		size_t u1 = m_cf[sep] + m_occ.rank(sep, u);
		if (l1 < u1)
			*out++ = FMInterval(l1, u1, it - first, last - first);
	}
	return out;
}

/** Search for a suffix of the query that matches a prefix of the
 * target.
 */
template <typename OutIt>
OutIt findOverlapSuffix(const std::string& q, OutIt out,
		unsigned minOverlap) const
{
	std::string s = q;
	std::transform(s.begin(), s.end(), s.begin(), Translate(*this));
	return findOverlapSuffix(s.begin(), s.end(), out, minOverlap);
}

/** Search for a prefix of the query that matches a suffix of the
 * target.
 */
template <typename It, typename OutIt>
OutIt findOverlapPrefix(It first, It last, OutIt out) const
{
	assert(first < last);

	size_t l = 1, u = m_occ.size();
	const char sep = 0;
	l = m_cf[sep] + m_occ.rank(sep, l);
	u = m_cf[sep] + m_occ.rank(sep, u);

	for (It it = last - 1; it >= first && l < u; --it) {
		T c = *it;
		if (c == SENTINEL())
			break;
		l = m_cf[c] + m_occ.rank(c, l);
		u = m_cf[c] + m_occ.rank(c, u);
	}
	if (l < u)
		*out++ = FMInterval(l, u, 0, last - first);
	return out;
}

/** Search for a prefix of the query that matches a suffix of the
 * target.
 */
template <typename OutIt>
OutIt findOverlapPrefix(const std::string& q, OutIt out,
		unsigned minOverlap) const
{
	std::string s = q;
	std::transform(s.begin(), s.end(), s.begin(), Translate(*this));
	typedef std::string::const_iterator It;
	It first = s.begin();
	for (It it = first + minOverlap; it <= s.end(); ++it)
		out = findOverlapPrefix(first, it, out);
	return out;
}

/** Search for a matching suffix of the query. */
template <typename It, typename MemoIt>
FMInterval findSuffix(It first, It last, MemoIt memoIt) const
{
	assert(first < last);
	typedef typename std::iterator_traits<MemoIt>::value_type
		Interval;

	size_t l = 1, u = m_occ.size();
	It it;
	for (it = last - 1; it >= first && l < u; --it) {
		T c = *it;
		if (c == SENTINEL())
			break;
		size_t l1 = m_cf[c] + m_occ.rank(c, l);
		size_t u1 = m_cf[c] + m_occ.rank(c, u);
		if (l1 >= u1)
			break;
		l = l1;
		u = u1;

		Interval fmi = Interval(l, u);
		if (*memoIt == fmi) {
			// This vertex of the prefix DAWG has been visited.
			break;
		}
		*memoIt++ = fmi;
	}
	return FMInterval(l, u, it - first + 1, last - first);
}

/** Search for a matching substring of the query at least k long.
 * @return the longest match
 */
template <typename It>
FMInterval findSubstring(It first, It last, unsigned k) const
{
	assert(first < last);
	FMInterval best(0, 0, 0, k > 0 ? k - 1 : 0);

	// Record which vertices of the prefix DAWG have been visited.
	std::vector<std::pair<size_type, size_type> > memo(last - first);
	std::pair<size_type, size_type>* memoIt = memo.data();
	for (It it = last; it > first; --it) {
		if (unsigned(it - first) < best.qspan())
			return best;
		FMInterval interval = findSuffix(first, it, memoIt++);
		if (interval.qspan() > best.qspan())
			best = interval;
	}
	return best;
}

/** Translate from ASCII to the indexed alphabet. */
struct Translate {
	Translate(const FMIndex& fmIndex) : m_fmIndex(fmIndex) { }
	T operator()(unsigned char c) const
	{
		return c < m_fmIndex.m_mapping.size()
			? m_fmIndex.m_mapping[c] : SENTINEL();
	}
  private:
	const FMIndex& m_fmIndex;
};

/** Search for a matching substring of the query at least k long.
 * @return the longest match and the number of matches
 */
Match find(const std::string& q, unsigned k) const
{
	std::string s = q;
	std::transform(s.begin(), s.end(), s.begin(), Translate(*this));

	FMInterval interval = findSubstring(s.begin(), s.end(), k);
	assert(interval.l <= interval.u);
	size_t count = interval.u - interval.l;
	if (count == 0)
		return Match(0, 0, 0, 0);
	return Match(interval.qstart, interval.qend,
			locate(interval.l), count);
}

/** Set the alphabet to [first, last). */
template <typename It>
void setAlphabet(It first, It last)
{
	std::vector<bool> mask(std::numeric_limits<T>::max());
	for (It it = first; it < last; ++it) {
		assert((size_t)*it < mask.size());
		mask[*it] = true;
	}

	m_alphabet.clear();
	m_mapping.clear();
	for (unsigned c = 0; c < mask.size(); ++c) {
		if (!mask[c])
			continue;
		m_mapping.resize(c + 1, std::numeric_limits<T>::max());
		m_mapping[c] = m_alphabet.size();
		m_alphabet.push_back(c);
	}
	assert(!m_alphabet.empty());
}

/** Set the alphabet to the characters of s. */
void setAlphabet(const std::string& s)
{
	setAlphabet(s.begin(), s.end());
}

#define STRINGIFY(X) #X
#define FM_VERSION_BITS(BITS) "FM " STRINGIFY(BITS) " 1"
#define FM_VERSION FM_VERSION_BITS(FMBITS)

/** Store an index. */
friend std::ostream& operator<<(std::ostream& out, const FMIndex& o)
{
	out << FM_VERSION << '\n'
		<< o.m_sampleSA << '\n';

	out << o.m_alphabet.size() << '\n';
	out.write((char*)o.m_alphabet.data(),
			o.m_alphabet.size() * sizeof o.m_alphabet[0]);

	out << o.m_sa.size() << '\n';
	out.write((char*)o.m_sa.data(),
		o.m_sa.size() * sizeof o.m_sa[0]);

	return out << o.m_occ;
}

/** Load an index. */
friend std::istream& operator>>(std::istream& in, FMIndex& o)
{
	std::string version;
	std::getline(in, version);
	assert(in);
	if (version != FM_VERSION) {
		std::cerr << "error: the version of this FM-index, `"
			<< version << "', does not match the version required "
			"by this program, `" FM_VERSION "'.\n";
		exit(EXIT_FAILURE);
	}

	in >> o.m_sampleSA;
	assert(in);

	size_t n;
	char c;

	in >> n;
	assert(in);
	c = in.get();
	assert(c == '\n');
	assert(n < std::numeric_limits<size_type>::max());
	o.m_alphabet.resize(n);
	in.read((char*)o.m_alphabet.data(), n * sizeof o.m_alphabet[0]);
	o.setAlphabet(o.m_alphabet.begin(), o.m_alphabet.end());

	in >> n;
	assert(in);
	c = in.get();
	assert(c == '\n');
	assert(n < std::numeric_limits<size_type>::max());
	o.m_sa.resize(n);
	in.read((char*)o.m_sa.data(), n * sizeof o.m_sa[0]);

	in >> o.m_occ;
	assert(in);
	o.countOccurences();

	return in;
}

private:

/** Build the cumulative frequency table m_cf from m_occ. */
void countOccurences()
{
	assert(!m_alphabet.empty());
	m_cf.resize(m_alphabet.size());
	// The sentinel character occurs once.
	m_cf[0] = 1;
	for (unsigned i = 0; i < m_cf.size() - 1; ++i)
		m_cf[i + 1] = m_cf[i] + m_occ.count(i);
}

	unsigned m_sampleSA;
	std::vector<T> m_alphabet;
	std::vector<T> m_mapping;
	std::vector<size_type> m_cf;
	std::vector<size_type> m_sa;
	BitArrays m_occ;
};

#endif
