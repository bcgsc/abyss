#include "Stats.h"
#include <math.h>

const double MINP = 0.00001f;

Stats::Stats() : m_maxDist(0)
{
	
}

//
// Perform a maximum likelihood estimate over the pdf and input distribution
//
int Stats::MaxLikelihoodEst(std::vector<int>& pairDistance)
{
	double maxL = -999999;
	int bestDist = -50;
	
	// TODO: Unhardcode this
	int minDist = -50;
	int maxDist = 200;
	for(int i = minDist; i < maxDist; i++)
	{
		double v = 	ComputeLikelihood(i, pairDistance);
		if(v > maxL)
		{
			maxL = v;
			bestDist = i;
		}
	}
	
	return bestDist;
}

//
// Compute the log likelihood function over the test distribution
//
double Stats::ComputeLikelihood(int d, std::vector<int>& testDist)
{
	double sum = 0.0f;

	for(std::vector<int>::iterator iter = testDist.begin(); iter != testDist.end(); iter++)
	{
		int val = *iter + d;
		double p;
		
		// Is this value in range of the pdf?
		if(val >= 0 && val < (int)m_pdf.size())
		{
			p = m_pdf[val];	
		}
		else
		{
			p = MINP;
		}
		
		sum += log(p);
	}
	
	return sum;
}

//
// Estimate number of pairs that would exist between two contigs given an (estimated distance) and the length of the second (non-root) contig
// The first contig is assumed to be long enough to not be limiting
//
int Stats::EstimateNumberOfPairs(int estDist, int contigLength, int numberOfReads)
{
	double p = 0.0f;
	for(int i = estDist; i < estDist + contigLength && i < m_maxDist; i++)
	{
		p += m_pdf[i];
	}
	
	assert(p <= 1.0f);
	return (int)round(p * (double)numberOfReads);
}

//
// Convert the histogram to a probability density function
// DOUBLE COUNTING PAIRS??
void Stats::GenerateStatsFromHistogram(const histogram& h)
{
	m_maxDist = 0;
	double count = 0;
	double sum = 0;
	for(histogram::const_iterator histIter = h.begin(); histIter != h.end(); histIter++)
	{
		if(histIter->first > m_maxDist)
		{
			m_maxDist = histIter->first;
		}
		count += histIter->second;
		sum += (histIter->second * histIter->first);
	}
	
	double mean = sum / count;
	
	// Compute the standard devition
	sum = 0;
	count = 0;
	for(histogram::const_iterator histIter = h.begin(); histIter != h.end(); histIter++)
	{

		count += histIter->second;
		sum += (static_cast<double>(histIter->second) * pow((static_cast<double>(histIter->first) - mean), 2));
	}
	
	sum /= count;
	m_stdDev = sqrt(sum);
	m_mean = mean;
	printf("mean distance: %lf std dev: %lf count: %d\n", m_mean, m_stdDev, (int)count);
	
	// Create the initial pdf with all values being 0
	m_pdf = PDF(m_maxDist+1, 0.0f);
	for(int i = 0; i < m_maxDist+1; i++)
	{
		histogram::const_iterator iter = h.find(i);
		if(iter != h.end())
		{
			m_pdf[i] = static_cast<double>(iter->second) / count;
		}
		else
		{
			m_pdf[i] = MINP;
		}
	}
}

//
// Get the standard deviation of the estimate based on the empirical standard deviation and the number of data points
//
double Stats::GetStdDevOfEstimate(int n)
{
	return m_stdDev / sqrt((double)n);
}
