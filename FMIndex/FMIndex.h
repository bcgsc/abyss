#ifndef FMINDEX_H
#define FMINDEX_H 1

#include "wat_array.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <istream>
#include <map>
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
};

/** An FM index. */
class FMIndex
{
  public:
	int load(std::istream& is);
	int save(std::ostream& os);
	int read(const char *fname, std::vector<uint8_t> &s);
	int buildFmIndex(const char *fnmae, int _percent);
	size_t locate(uint64_t i) const;

	FMIndex() : percent(0), alphaSize(0) { }

	/** Construct an FMIndex. */
	FMIndex(const std::string& path) : percent(0), alphaSize(0)
	{
		std::ifstream in(path.c_str());
		assert(in);
		int status = load(in);
		assert(status == 0);
	}

	/** Search for a matching suffix of the query. */
	template <typename It>
	std::pair<size_t, size_t> findExact(It first, It last) const
	{
		assert(first < last);
		size_t l = 1, u = wa.length();
		It it;
		for (it = last - 1; it >= first && l < u; --it) {
			uint8_t c = *it;
			l = cf[c] + wa.Rank(c, l);
			u = cf[c] + wa.Rank(c, u);
		}
		return std::make_pair(l, u);
	}

	/** Search for a matching suffix of the query. */
	template <typename It>
	FMInterval findSuffix(It first, It last) const
	{
		assert(first < last);
		size_t l = 1, u = wa.length();
		It it;
		for (it = last - 1; it >= first && l < u; --it) {
			uint8_t c = *it;
			size_t l1 = cf[c] + wa.Rank(c, l);
			size_t u1 = cf[c] + wa.Rank(c, u);
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
			if (unsigned(it - first) < best.qend - best.qstart)
				return best;
			FMInterval interval = findSuffix(first, it);
			if (interval.qend - interval.qstart
					> best.qend - best.qstart)
				best = interval;
		}
		return best;
	}

	/** Translate from ASCII to the indexed alphabet. */
	struct Translate {
		Translate(const FMIndex& fmIndex) : m_fmIndex(fmIndex) { }
		uint8_t operator()(unsigned char c) const
		{
			std::map<unsigned char, uint8_t>::const_iterator
				it = m_fmIndex.mapping.find(c);
			return it == m_fmIndex.mapping.end() ? m_fmIndex.alphaSize
				: it->second;
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

  private:
	void calculateStatistics(const std::vector<uint8_t> &s);
	int buildBWT(const std::vector<uint8_t> &s,
			const std::vector<uint32_t> &sa,
			std::vector<uint64_t> &bwt);
	int buildSA(const std::vector<uint8_t> &s,
			std::vector<uint32_t> &sa);
	int buildSampledSA(const std::vector<uint8_t> &s,
			const std::vector<uint32_t> &sa);
	int buildWaveletTree(const std::vector<uint64_t> &bwt);

	int percent;
	uint8_t alphaSize;
	std::vector<uint32_t> cf;
	std::map<unsigned char, uint8_t> mapping;
	std::map<uint8_t, unsigned char> rmapping;
	std::vector<uint32_t> sampledSA;
	std::vector<std::vector<int> > cols;
	wat_array::WatArray wa;
};

#endif
