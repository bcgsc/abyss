#ifndef STATS_H
#define STATS_H 1

#include <cmath>
#include <map>
#include <ostream>
#include <vector>

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

typedef std::vector<double> DoubleVec;
struct PDF
{
	PDF() {};
	PDF(const Histogram& h);
	
	double getP(size_t idx) const;
	double getMinP() const { return m_minp; }
	size_t getMaxIdx() const { return m_maxIdx; }
	void print() const;

	double getSampleStdDev(unsigned n) const
	{
		return m_stdDev / sqrt((double)n);
	}

	size_t m_maxIdx;
	DoubleVec m_dist;
	double m_mean;
	double m_stdDev;
	double m_minp;
	
	// calculate the minimal range in which p% of the values will fall into
	void calculateMinimalRange(double p, size_t& low, size_t& high) const;
};

// Maximum Likelihood Estimator functions
int maxLikelihoodEst(int min, int max,
		const std::vector<int>& pairDistance, const PDF& pdf,
		unsigned& n);

// Compute the likelihood of the distribution
double computeLikelihood(int param, const std::vector<int>& testDist,
		const PDF& pdf, unsigned& n);

#endif
