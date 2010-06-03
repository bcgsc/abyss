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

// Maximum Likelihood Estimator functions
int maxLikelihoodEst(int min, int max,
		const std::vector<int>& pairDistance, const PDF& pdf,
		unsigned& n);

// Compute the likelihood of the distribution
double computeLikelihood(int param, const std::vector<int>& testDist,
		const PDF& pdf, unsigned& n);

#endif
