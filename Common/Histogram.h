#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1

#include <map>
#include <vector>
#include <ostream>

/** A histogram of type T, which is int be default.
 * A histogram may be implemented as a multiset. This class aims
 * to provide a similar interface to a multiset.
 */
class Histogram : public std::map<int, unsigned>
{
  public:
	typedef int T;
	typedef std::map<T, unsigned> Map;

	Histogram() {}
	Histogram(const std::vector<T>& data)
	{
		for (std::vector<T>::const_iterator iter = data.begin();
				iter != data.end(); ++iter)
			insert(*iter);
	}

	void insert(T value) { (*this)[value]++; }

	void insert(T value, unsigned count)
	{
		(*this)[value] += count;
	}

	unsigned size() const
	{
		T min = 0;
		T max = maximum();
		unsigned sum = 0;
		for (T value = min; value <= max; ++value)
			sum += count(value);
		return sum;
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
