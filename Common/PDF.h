#ifndef STATS_H
#define STATS_H 1

#include "Histogram.h"
#include <cassert>
#include <cmath>
#include <vector>

class Histogram;

struct PDF
{
	/** Construct a PDF from a histogram. */
	PDF(const Histogram& h) : m_dist(h.maximum() + 1), m_stdDev(h.sd())
	{
		unsigned count = h.size();
		m_minp = (double)1 / count;
		for (size_t i = 0; i < m_dist.size(); i++) {
			unsigned n = h.count(i);
			m_dist[i] = n > 0 ? (double)n / count : m_minp;
		}
	}

	/** Return the probability of x. */
	double getP(size_t x) const
	{
		return x < m_dist.size() ? m_dist[x] : m_minp;
	}

	double getMinP() const { return m_minp; }

	size_t getMaxIdx() const {
		assert(!m_dist.empty());
		return m_dist.size() - 1;
	}

	double getSampleStdDev(unsigned n) const
	{
		return m_stdDev / sqrt(n);
	}

	std::vector<double> m_dist;
	double m_stdDev;
	double m_minp;
};

#endif
