#ifndef PMF_H
#define PMF_H 1

#include "Common/Exception.h"
#include "Histogram.h"
#include <cassert>
#include <cmath>
#include <vector>

class Histogram;

/** Probability mass function */
class PMF
{
  public:
	/** Construct a PMF from a histogram. */
	PMF(const Histogram& h)
		: m_dist(h.maximum() + 1), m_mean(h.mean()), m_stdDev(h.sd()), m_median(h.median())
	{
		unsigned count = h.size();
		m_minp = (double)1 / count;
		for (size_t i = 0; i < m_dist.size(); i++) {
			unsigned n = h.count(i);
			m_dist[i] = n > 0 ? (double)n / count : m_minp;
		}
	}

	/** Return the probability of x. */
	double operator[](size_t x) const
	{
		return x < m_dist.size() ? m_dist[x] : m_minp;
	}

	/** Return the minimum probability. */
	double minProbability() const { return m_minp; }

	/** Return the minimum value. */
	size_t minValue() const { return 0; }

	/** Return the maximum value. */
	size_t maxValue() const
	{
		assert(!m_dist.empty());
		return m_dist.size() - 1;
	}

	/** Return the median of this distribution. */
	int median() const { return m_median; }

	/** Return the mean of this distribution. */
	double mean() const { return m_mean; }

	/** Return the standard deviation of the sampled mean
	 * of n observations.
	 */
	double getSampleStdDev(unsigned n) const
	{
		return m_stdDev / sqrt(n);
	}

  private:
	std::vector<double> m_dist;
	double m_mean;
	double m_stdDev;
	double m_minp;
	int m_median;
};

namespace std {
	template<>
	inline void swap(PMF&, PMF&) NOEXCEPT { assert(false); }
}

#endif
