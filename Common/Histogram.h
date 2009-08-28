#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1

#include <map>
#include <vector>
#include <ostream>

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
			addDataPoint(*iter);
	}

	void addDataPoint(T value) { (*this)[value]++; }

	void addMultiplePoints(T value, unsigned count)
	{
		(*this)[value] += count;
	}

	Histogram trim(double percent);

	unsigned getSumCount() const
	{
		T min = 0;
		T max = getMax();
		unsigned sum = 0;
		for (T value = min; value <= max; ++value)
			sum += getCount(value);
		return sum;
	}

	unsigned getCount(T value) const
	{
		Map::const_iterator iter = find(value);
		return find(value) == end() ? 0 : iter->second;
	}

	T getMin() const
	{
		return empty() ? 0 : begin()->first;
	}

	T getMax() const
	{
		return empty() ? 0 : rbegin()->first;
	}

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
