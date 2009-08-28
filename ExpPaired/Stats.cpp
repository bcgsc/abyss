#include "Stats.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

//
// Perform a maximum likelihood estimate over the pdf and input distribution
//
int maxLikelihoodEst(int min, int max,
		const std::vector<int>& pairDistance, const PDF& pdf,
		unsigned& n)
{
	double maxL = -999999;
	double nextBestL = maxL;
	int bestDist = min;
	unsigned bestn = 0;
	
	int minDist = min;
	int maxDist = max;
	for(int i = minDist; i < maxDist; i++)
	{
		unsigned trialn;
		double v = computeLikelihood(i, pairDistance, pdf, trialn);
		if(v > maxL)
		{
			nextBestL = maxL;
			maxL = v;
			bestDist = i;
			bestn = trialn;
		}
	}

	n = bestn;
	return bestDist;
}

//
// Compute the log likelihood function over the test distribution
//
double computeLikelihood(int param, const std::vector<int>& testDist,
		const PDF& pdf, unsigned &n)
{
	n = 0;
	double sum = 0.0f;
	for (std::vector<int>::const_iterator iter = testDist.begin();
			iter != testDist.end(); ++iter) {
		int val = *iter + param;
		double p = pdf.getP(val);
		sum += log(p);
		if (p > pdf.getMinP())
			n++;
	}
	return sum;
}

/** Trim off the bottom percent/2 and top percent/2 data points.
 * At least (1 - percent) of the data will remain.
 */
Histogram Histogram::trim(double percent)
{
	// The amount to take off each end
	double half_percent = percent/2;
	double low_cutoff = half_percent;
	double high_cutoff = 1.0f - half_percent;
	double total = (double)getSumCount();

	T min = 0;
	T max = getMax();
	double cumulative = 0.0f;

	Histogram newHist;
	for (T value = min; value <= max; ++value) {
		unsigned currCount = getCount(value);
		double frac = currCount / total;
		double temp_total = cumulative + frac;
		//printf("Frac: %lf TT: %lf C: %lf LC: %lf HC: %lf\n", frac, temp_total, cumulative, low_cutoff, high_cutoff);
		if(temp_total > low_cutoff && cumulative < high_cutoff)
		{
			// Add these elements
			newHist.addMultiplePoints(value, currCount);
		}
		cumulative = temp_total;
	}

	return newHist;
}

// Construct a pdf from a histogram
PDF::PDF(const Histogram& h)
{
	m_maxIdx = h.getMax();
	unsigned count = 0;

	for(Histogram::Map::const_iterator histIter = h.begin();
			histIter != h.end(); histIter++)
		count += histIter->second;
	m_minp = (double)1 / count;

	// Create the initial pdf with all values being 0
	m_dist = DoubleVec(m_maxIdx+1, 0.0f);
	for(size_t i = 0; i <= m_maxIdx; i++)
	{
		unsigned v = h.getCount(i);
		m_dist[i] = v > 0 ? (double)v / count : m_minp;
	}

	// Calculate the mean
	double sum = 0;
	double sumsqr = 0;
	for(Histogram::Map::const_iterator histIter = h.begin();
			histIter != h.end(); histIter++) {
		sum += ((double)histIter->second * (double)histIter->first);
		sumsqr += histIter->second * (pow(histIter->first, 2.0f));
	}

	m_mean = sum / count;
	assert(m_mean > 0);
	double msqr = pow(sum,2.0) / count;
	double ss1 = sumsqr - msqr;
	double t1 = ss1 / count;

	m_stdDev = sqrt(t1);
	
	printf("Stats mean: %.2lf sd: %.2lf n: %u min: %u max: %u\n",
			m_mean, m_stdDev, count, h.getMin(), h.getMax());
}

double PDF::getP(size_t idx) const
{
	return (idx <= m_maxIdx) ? m_dist[idx] : m_minp;
}

//
//
//
void PDF::print() const
{
	for(size_t idx = 0; idx <= m_maxIdx; ++idx)
	{
		printf("%zu %lf\n", idx, m_dist[idx]);
	}
}

//
//
//
void PDF::calculateMinimalRange(double p, size_t& low, size_t& high) const
{
	// First, find the max prob
	double max = 0.0f;
	size_t maxIdx = 0;
	for(size_t idx = 0; idx <= m_maxIdx; ++idx)
	{
		if(m_dist[idx] > max)
		{
			max = m_dist[idx];
			maxIdx = idx;
		}
	}
	
	size_t leftIdx = maxIdx;
	size_t rightIdx = maxIdx;
	
	low = 0;
	high = 0;
	double cum = m_dist[maxIdx];
	while(cum < p)
	{		
		// compute the next index
		size_t nextLeft = (leftIdx > 0) ? leftIdx - 1 : leftIdx;
		size_t nextRight = (rightIdx < m_maxIdx) ? rightIdx + 1 : rightIdx;
		
		// If the range could not be expanded, terminate
		if(nextLeft == leftIdx && nextRight == rightIdx)
		{
			break;
		}
		else
		{
			// left expansion gives the greatest improvement
			if((m_dist[nextLeft] > m_dist[nextRight] && nextLeft != leftIdx) || (nextRight == rightIdx))
			{
				cum += m_dist[nextLeft];
				leftIdx = nextLeft;
			}
			else // right expansion gives the greatest improvement
			{
				cum += m_dist[nextRight];
				rightIdx = nextRight;
			}
		}
	}
	
	low = leftIdx;
	high = rightIdx;
}
