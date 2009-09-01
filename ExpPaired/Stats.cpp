#include "Stats.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

// Perform a maximum likelihood estimate over the pdf and input distribution
int maxLikelihoodEst(int min, int max,
		const vector<int>& pairDistance, const PDF& pdf, unsigned& n)
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
double computeLikelihood(int param, const vector<int>& testDist,
		const PDF& pdf, unsigned &n)
{
	n = 0;
	double sum = 0.0f;
	for (vector<int>::const_iterator iter = testDist.begin();
			iter != testDist.end(); ++iter) {
		int val = *iter + param;
		double p = pdf.getP(val);
		sum += log(p);
		if (p > pdf.getMinP())
			n++;
	}
	return sum;
}

// Construct a pdf from a histogram
PDF::PDF(const Histogram& h) :
	m_maxIdx(h.maximum()),
	m_dist(m_maxIdx + 1),
	m_mean(h.mean()),
	m_stdDev(h.sd())
{
	unsigned count = h.size();
	m_minp = (double)1 / count;

	for (size_t i = 0; i <= m_maxIdx; i++) {
		unsigned v = h.count(i);
		m_dist[i] = v > 0 ? (double)v / count : m_minp;
	}

	printf("Stats mean: %.2lf sd: %.2lf n: %u min: %u max: %u\n",
			m_mean, m_stdDev, count, h.minimum(), h.maximum());
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
