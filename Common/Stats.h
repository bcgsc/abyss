#ifndef STATS_H
#define STATS_H 1

#include <cmath>
#include <vector>

class Histogram;

struct PDF
{
	PDF() {};
	PDF(const Histogram& h);
	
	double getP(size_t idx) const;
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

int maximumLikelihoodEstimate(int first, int last,
		const std::vector<int>& samples, const PDF& pdf,
		unsigned len0, unsigned len1, unsigned& n);

#endif
