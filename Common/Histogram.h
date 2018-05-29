#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1

#include "Common/Exception.h"
#include "StringUtil.h" // for toEng
#include "VectorUtil.h" // for make_vector
#include <cassert>
#include <climits> // for INT_MAX
#include <cmath>
#include <istream>
#include <map>
#include <ostream>
#include <string>
#include <vector>

/** A histogram of type T, which is int be default.
 * A histogram may be implemented as a multiset. This class aims
 * to provide a similar interface to a multiset.
 */
class Histogram
{
	typedef int T;
	typedef size_t size_type;
	typedef std::map<T, size_type> Map;
	typedef long long unsigned accumulator;

  public:
	typedef Map::const_iterator const_iterator;

	Histogram() { }

	/** Construct a histogram of the specified elements. */
	template <class InputIterator>
	Histogram(InputIterator first, InputIterator last)
	{
		for (InputIterator it = first; it != last; ++it)
			insert(*it);
	}

	/** Construct a histogram from a vector, where the index into the
	 * vector is the sample, and the value at that index is the number
	 * of times that sample was observed.
	 */
	explicit Histogram(std::vector<size_type> v)
	{
		for (T i = 0; i < (T)v.size(); i++)
			if (v[i] > 0)
				m_map.insert(std::make_pair(i, v[i]));
	}

	void insert(T value) { m_map[value]++; }

	void insert(T value, size_type count) { m_map[value] += count; }

	size_type count(T value) const
	{
		Map::const_iterator iter = m_map.find(value);
		return iter == m_map.end() ? 0 : iter->second;
	}

	/** Return the number of elements in the range [lo,hi]. */
	size_type count(T lo, T hi) const
	{
		assert(lo <= hi);
		size_type n = 0;
		Map::const_iterator last = m_map.upper_bound(hi);
		for (Map::const_iterator it = m_map.lower_bound(lo);
				it != last; ++it)
			n += it->second;
		return n;
	}

	T minimum() const
	{
		return empty() ? 0 : m_map.begin()->first;
	}

	T maximum() const
	{
		return empty() ? 0 : m_map.rbegin()->first;
	}

	bool empty() const { return m_map.empty(); }

	size_type size() const
	{
		size_type n = 0;
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it)
			n += it->second;
		return n;
	}

	/** Return the sum. */
	accumulator sum() const
	{
		accumulator total = 0;
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it)
			total += (accumulator)it->first * it->second;
		return total;
	}

	/** Return the mean. */
	double mean() const
	{
		accumulator n = 0, total = 0;
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it) {
			n += it->second;
			total += (accumulator)it->first * it->second;
		}
		return (double)total / n;
	}

	double variance() const
	{
		accumulator n = 0, total = 0, squares = 0;
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it) {
			n += it->second;
			total += (accumulator)it->first * it->second;
			squares += (accumulator)it->first * it->first
				* it->second;
		}
		return (squares - (double)total * total / n) / n;
	}

	double sd() const
	{
		return sqrt(variance());
	}

	/** Return the specified percentile. */
	T percentile(float p) const
	{
		size_type x = (size_type)ceil(p * size());
		size_type n = 0;
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it) {
			n += it->second;
			if (n >= x)
				return it->first;
		}
		return maximum();
	}

	/** Return the median. */
	T median() const
	{
		return percentile(0.5);
	}

	/** Return the largest weight in the arg min of partial sum of
	 * weights. */
	T argMin(accumulator x) const
	{
		accumulator total = 0;
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it) {
			total += (accumulator)it->first * it->second;
			if (total >= x)
				return it->first;
		}
		return maximum();
	}

	/** Return the specified weighted percentile. */
	T weightedPercentile(float p) const
	{
		return argMin((accumulator)ceil(p * sum()));
	}

	/** Return the expected value */
	double expectedValue() const
	{
		double value = 0;
		accumulator acc = sum();
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); it++) {
			value += (double)it->first * it->first
				* it->second / acc;
		}
		return value;
	}

	/** Return the N50. */
	T n50() const { return weightedPercentile(0.5); }

	/** Return the first local minimum or zero if a minimum is not
	 * found. */
	T firstLocalMinimum() const
	{
		const unsigned SMOOTHING = 4;
		assert(!empty());
		Map::const_iterator minimum = m_map.begin();
		size_type count = 0;
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it) {
			if (it->second <= minimum->second) {
				minimum = it;
				count = 0;
			} else if (++count >= SMOOTHING)
				break;
		}
		if (minimum->first == maximum())
			return 0;
		return minimum->first;
	}

	void eraseNegative()
	{
		for (Map::iterator it = m_map.begin(); it != m_map.end();)
			if (it->first < 0)
				m_map.erase(it++);
			else
				++it;
	}

	/** Remove noise from the histogram. Noise is defined as a
	 * sample x where h[x-1] == 0 && h[x+1] == 0.
	 */
	void removeNoise()
	{
		for (Map::iterator it = m_map.begin(); it != m_map.end();) {
			if (m_map.count(it->first - 1) == 0
					&& m_map.count(it->first + 1) == 0
					&& m_map.size() > 1)
				m_map.erase(it++);
			else
				++it;
		}
	}

	/** Remove outliers from the histogram. A sample is an outlier
	 * if it is outside the range [Q1 - k*(Q3-Q1), Q3 + k*(Q3-Q1)]
	 * where k = 20.
	 */
	void removeOutliers()
	{
		T q1 = percentile(0.25);
		T q3 = percentile(0.75);
		T l = q1 - 20 * (q3 - q1);
		T u = q3 + 20 * (q3 - q1);
		for (Map::iterator it = m_map.begin(); it != m_map.end();) {
			if (it->first < l || it->first > u)
				m_map.erase(it++);
			else
				++it;
		}
	}

	/** Negate each element of this histogram. */
	Histogram negate() const
	{
		Histogram h;
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it)
			h.m_map.insert(std::make_pair(-it->first, it->second));
		return h;
	}

	Histogram trimFraction(double fraction) const;
	Histogram trimLow(T threshold) const;

	typedef std::vector<accumulator> Bins;
	Bins bin(unsigned n) const;
	std::string barplot() const;
	std::string barplot(unsigned nbins) const;

	const_iterator begin() const { return m_map.begin(); }
	const_iterator end() const { return m_map.end(); }

	/** Return a vector representing this histogram. */
	std::vector<size_type> toVector() const
	{
		assert(minimum() >= 0);
#if 0
		std::vector<size_type> v(maximum()+1);
#else
		// CommLayer::reduce requires the arrays have the same size.
		std::vector<size_type> v(65536);
		assert(maximum() < (T)v.size());
#endif
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it)
			v[it->first] = it->second;
		return v;
	}

	friend std::ostream& operator<<(std::ostream& o,
			const Histogram& h)
	{
		for (Map::const_iterator it = h.m_map.begin();
				it != h.m_map.end(); ++it)
			o << it->first << '\t' << it->second << '\n';
		return o;
	}

	friend std::istream& operator>>(std::istream& in, Histogram& h)
	{
		Histogram::T value;
		size_type count;
		while (in >> value >> count)
			h.insert(value, count);
		assert(in.eof());
		return in;
	}

  private:
	Map m_map;
};

namespace std {
	template<>
	inline void swap(Histogram&, Histogram&) NOEXCEPT { assert(false); }
}

/** Print assembly contiguity statistics header. */
static inline std::ostream& printContiguityStatsHeader(
		std::ostream& out,
		unsigned minSize,
		const std::string& sep = "\t",
		const long long unsigned expSize = 0)
{
	out << "n" << sep
		<< "n:" << minSize << sep
		<< "L50" << sep;
	if (expSize > 0)
		out << "LG50" << sep
			<< "NG50" << sep;
	return out << "min" << sep
		<< "N75" << sep
		<< "N50" << sep
		<< "N25" << sep
		<< "E-size" << sep
		<< "max" << sep
		<< "sum" << sep
		<< "name" << '\n';
}

/** Print assembly contiguity statistics. */
static inline std::ostream& printContiguityStats(
		std::ostream& out, const Histogram& h0,
		unsigned minSize, bool printHeader = true,
		const std::string& sep = "\t",
		const long long unsigned expSize = 0)
{
	Histogram h = h0.trimLow(minSize);
	if (printHeader)
		printContiguityStatsHeader(out, minSize, sep, expSize);
	unsigned n50 = h.n50();
	out << toEng(h0.size()) << sep
		<< toEng(h.size()) << sep
		<< toEng(h.count(n50, INT_MAX)) << sep;
	long long unsigned sum = h.sum();
	if (expSize > 0) {
		unsigned ng50;
		if (sum < expSize/2)
			ng50 = h.minimum();
		else
			ng50 = h.argMin(sum - expSize/2);
		out << toEng(h.count(ng50, INT_MAX)) << sep
			<< toEng(ng50) << sep;
	}
	return out
		<< toEng(h.minimum()) << sep
		<< toEng(h.weightedPercentile(1 - 0.75)) << sep
		<< toEng(n50) << sep
		<< toEng(h.weightedPercentile(1 - 0.25)) << sep
		<< toEng((unsigned)h.expectedValue()) << sep
		<< toEng(h.maximum()) << sep
		<< toEng(sum);
}

/** Pass assembly contiguity statistics -- values only. */
static inline std::vector<int> passContiguityStatsVal(
		const Histogram& h0, unsigned minSize, const long long unsigned expSize = 0)
{
#if _SQL
	Histogram h = h0.trimLow(minSize);
	unsigned n50 = h.n50();
	long long unsigned sum = h.sum();

	std::vector<int> vec = make_vector<int>()
		<< h0.size()
		<< h.size()
		<< h.count(n50, INT_MAX)
		<< h.minimum()
		<< h.weightedPercentile(1 - 0.75)
		<< n50
		<< h.weightedPercentile(1 - 0.25)
		<< (unsigned)h.expectedValue()
		<< h.maximum()
		<< sum;

	if (expSize > 0) {
		unsigned ng50;
		if (sum < expSize/2)
			ng50 = h.minimum();
		else
			ng50 = h.argMin(sum - expSize/2);
		vec.push_back(h.count(ng50, INT_MAX));
		vec.push_back(ng50);
	}

	return vec;
#else
	(void)h0;
	(void)minSize;
	(void)expSize;
	return make_vector<int>();
#endif
}
#endif
