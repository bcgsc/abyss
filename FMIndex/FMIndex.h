#ifndef FMINDEX_H
#define FMINDEX_H 1

#include "BitArrays.h"
#include <algorithm>
#include <cassert>
#include <climits> // for UCHAR_MAX
#include <fstream>
#include <istream>
#include <limits> // for numeric_limits
#include <ostream>
#include <stdint.h>
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
	/** An index. */
	typedef uint32_t size_type;

	/** A symbol. */
	typedef uint8_t T;

	/** The sentinel symbol. */
	static T SENTINEL() { return std::numeric_limits<T>::max(); }

  public:
	void load(std::istream& is);
	void save(std::ostream& os);
	void read(const char *path, std::vector<T> &s);
	void buildFmIndex(const char *path, unsigned sampleSA);
	size_t locate(size_t i) const;

	FMIndex() : m_sampleSA(0), m_alphaSize(0) { }

	/** Construct an FMIndex. */
	FMIndex(const std::string& path) : m_sampleSA(0), m_alphaSize(0)
	{
		std::ifstream in(path.c_str());
		assert(in);
		load(in);
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

		// Construct the reverse transformation of the alphabet.
		std::vector<T> decode(m_alphaSize, SENTINEL());
		for (unsigned i = 0; i < m_mapping.size(); ++i) {
			T c = m_mapping[i];
			if (c == SENTINEL())
				continue;
			assert(c < decode.size());
			decode[c] = i;
		}

		// Translate the character set and output the result.
		for (std::vector<T>::const_reverse_iterator it = s.rbegin();
				it != s.rend(); ++it) {
			T c = *it;
			assert(c < decode.size());
			*out++ = decode[c];
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
			if (c == UCHAR_MAX)
				return std::make_pair(0, 0);
			l = m_cf[c] + m_occ.rank(c, l);
			u = m_cf[c] + m_occ.rank(c, u);
		}
		return std::make_pair(l, u);
	}

	/** Search for an exact match. */
	template <typename T>
	std::pair<size_t, size_t> findExact(const T& q) const
	{
		T s(q.size());
		std::transform(q.begin(), q.end(), s.begin(),
				Translate(*this));
		return findExact(s.begin(), s.end());
	}

	/** Search for a matching suffix of the query. */
	template <typename It>
	FMInterval findSuffix(It first, It last) const
	{
		assert(first < last);
		size_t l = 1, u = m_occ.size();
		It it;
		for (it = last - 1; it >= first && l < u; --it) {
			T c = *it;
			if (c == UCHAR_MAX)
				break;
			size_t l1 = m_cf[c] + m_occ.rank(c, l);
			size_t u1 = m_cf[c] + m_occ.rank(c, u);
			if (l1 >= u1)
				break;
			l = l1;
			u = u1;
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
		for (It it = last; it > first; --it) {
			if (unsigned(it - first) < best.qspan())
				return best;
			FMInterval interval = findSuffix(first, it);
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
				? m_fmIndex.m_mapping[c] : UCHAR_MAX;
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
		std::transform(s.begin(), s.end(), s.begin(),
				Translate(*this));

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
		assert(m_mapping.empty() && m_alphaSize == 0);
		for (It it = first; it < last; ++it) {
			unsigned c = *it;
			if (c >= m_mapping.size())
				m_mapping.resize(c + 1, UCHAR_MAX);
			if (m_mapping[c] == UCHAR_MAX) {
				m_mapping[c] = m_alphaSize++;
				assert(m_alphaSize < UCHAR_MAX);
			}
		}
		assert(m_alphaSize > 0);
	}

	/** Set the alphabet to the characters of s. */
	void setAlphabet(const std::string& s)
	{
		setAlphabet(s.begin(), s.end());
	}

  private:
	void calculateStatistics(const std::vector<T> &s);
	void buildBWT(const std::vector<T> &s,
			const std::vector<size_type> &sa,
			std::vector<T> &bwt);
	void buildSA(const std::vector<T> &s,
			std::vector<size_type> &sa);
	void buildSampledSA(const std::vector<size_type> &sa);

	unsigned m_sampleSA;
	unsigned m_alphaSize;
	std::vector<size_type> m_cf;
	std::vector<T> m_mapping;
	std::vector<size_type> m_sampledSA;
	BitArrays m_occ;
};

#endif
