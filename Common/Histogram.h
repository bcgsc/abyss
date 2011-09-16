#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1

#include <cassert>
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

	T median() const
	{
		size_type half = (size() + 1) / 2;
		size_type n = 0;
		for (Map::const_iterator it = m_map.begin();
				it != m_map.end(); ++it) {
			n += it->second;
			if (n >= half)
				return it->first;
		}
		return maximum();
	}

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
	inline void swap(Histogram&, Histogram&) { assert(false); }
}

#endif
