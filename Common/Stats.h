#ifndef STATS_H
#define STATS_H 1

#include "Histogram.h"
#include <cmath>
#include <vector>

class Histogram;

struct PDF
{
	/** Construct a PDF from a histogram. */
	PDF(const Histogram& h) : m_maxIdx(h.maximum()), m_dist(m_maxIdx + 1),
		m_mean(h.mean()), m_stdDev(h.sd())
	{
		unsigned count = h.size();
		m_minp = (double)1 / count;

		for (size_t i = 0; i <= m_maxIdx; i++) {
			unsigned v = h.count(i);
			m_dist[i] = v > 0 ? (double)v / count : m_minp;
		}
	}

	/** Return the probability of x. */
	double getP(size_t x) const
	{
		return x <= m_maxIdx ? m_dist[x] : m_minp;
	}

	double getMinP() const { return m_minp; }
	size_t getMaxIdx() const { return m_maxIdx; }

	double getSampleStdDev(unsigned n) const
	{
		return m_stdDev / sqrt((double)n);
	}

	size_t m_maxIdx;
	std::vector<double> m_dist;
	double m_mean;
	double m_stdDev;
	double m_minp;

	// calculate the minimal range in which p% of the values will fall into
	void calculateMinimalRange(double p, size_t& low, size_t& high) const;
};

#endif
