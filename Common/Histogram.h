#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1

#include <cmath>
#include <map>
#include <ostream>

/** A histogram of type T, which is int be default.
 * A histogram may be implemented as a multiset. This class aims
 * to provide a similar interface to a multiset.
 */
class Histogram : std::map<int, unsigned>
{
  public:
	typedef int T;
	typedef std::map<T, unsigned> Map;

	void insert(T value) { (*this)[value]++; }

	void insert(T value, unsigned count)
	{
		(*this)[value] += count;
	}

	unsigned count(T value) const
	{
		Map::const_iterator iter = find(value);
		return find(value) == end() ? 0 : iter->second;
	}

	T minimum() const
	{
		return empty() ? 0 : begin()->first;
	}

	T maximum() const
	{
		return empty() ? 0 : rbegin()->first;
	}

	unsigned size() const
	{
		unsigned n = 0;
		for (Histogram::Map::const_iterator it = begin();
				it != end(); it++)
			n += it->second;
		return n;
	}

	double mean() const
	{
		unsigned n = 0, total = 0;
		for (Histogram::Map::const_iterator it = begin();
				it != end(); it++) {
			n += it->second;
			total += it->first * it->second;
		}
		return (double)total / n;
	}

	double variance() const
	{
		unsigned n = 0, total = 0, squares = 0;
		for (Histogram::Map::const_iterator it = begin();
				it != end(); it++) {
			n += it->second;
			total += it->first * it->second;
			squares += it->first * it->first * it->second;
		}
		return (squares - total * total / (double)n) / n;
	}

	double sd() const
	{
		return sqrt(variance());
	}

	void eraseNegative()
	{
		for (Histogram::Map::iterator it = begin();
				it != end();)
			if (it->first < 0)
				erase(it++);
			else
				++it;
	}

	Histogram trim(double percent) const;

	friend std::ostream& operator<<(std::ostream& o,
			const Histogram& h)
	{
		for (Map::const_iterator it = h.begin();
				it != h.end(); ++it)
			o << it->first << '\t' << it->second << '\n';
		return o;
	}
};

#endif
