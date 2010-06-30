#include "Stats.h"
#include "Histogram.h"
#include <algorithm>
#include <cassert>
#include <climits> // for INT_MAX
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

/** This window function is a triangle with a flat top, or a rectangle
 * with sloped sides.
 * @param len0 a parameter of the window function
 * @param len1 a parameter of the window function
 * @param d a parameter of the window function
 * @param x the value at which to evaluate the window function
 * @return the window function evaluated at x
 */
static double window(int len0, int len1, int d, int x)
{
	assert(len0 > 0);
	assert(len1 > 0);
	if (d < 0)
		d = 0;
	return (x <= d
			? 1
		: x < d + len0
			? x - d
		: x < d + len1
			? len0
		: x < d + len0 + len1
			? d + len0 + len1 - x
		: 1) / (double)len0;
}

/** Compute the log likelihood that these samples came from the
 * specified distribution scaled by a window function.
 * @param theta the parameter of the PDF
 * @param samples the samples
 * @param pdf the PDF scaled by the window function
 * @param len0 a parameter of the window function
 * @param len1 a parameter of the window function
 * @param d a parameter of the window function
 * @param n [out] the number of samples with a non-zero probability
 * @return the likelihood
 */
static double computeLikelihood(int theta,
		const vector<int>& samples, const PDF& pdf,
		unsigned len0, unsigned len1, int d,
		unsigned &n)
{
	assert(len0 <= len1);
	n = 0;
	double sum = 0.0f;
	for (vector<int>::const_iterator iter = samples.begin();
			iter != samples.end(); ++iter) {
		int val = *iter + theta;
		double p = pdf.getP(val);
		/* When randomly selecting fragments that span a given point,
		 * longer fragments are more likely to be selected than
		 * shorter fragments.
		 */
		sum += log(p * window(len0, len1, d, val));
		if (p > pdf.getMinP())
			n++;
	}
	return sum;
}

static int maximumLikelihoodEstimate(int first, int last,
		const vector<int>& samples, const PDF& pdf,
		unsigned len0, unsigned len1, int d,
		unsigned& n)
{
	double maxL = -999999;
	double nextBestL = maxL;
	int bestDist = first;
	unsigned bestn = 0;
	for (int theta = first; theta < last; theta++) {
		unsigned trialn;
		double v = computeLikelihood(theta, samples, pdf,
				len0, len1, d, trialn);
		if(v > maxL)
		{
			nextBestL = maxL;
			maxL = v;
			bestDist = theta;
			bestn = trialn;
		}
	}

	n = bestn;
	return bestDist;
}

/** Return the most likely distance between two contigs.
 * @param len0 the length of the first contig
 * @param len1 the length of the second contig
 */
int maximumLikelihoodEstimate(int first, int last,
		const vector<int>& samples, const PDF& pdf,
		unsigned len0, unsigned len1,
		unsigned& n)
{
	(void)len0; (void)len1;
	return maximumLikelihoodEstimate(first, last, samples, pdf,
				1, INT_MAX/2, 0, n);
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

	cerr << "Stats mean: " << setprecision(4) << m_mean << " "
	   "sd: " << setprecision(4) << m_stdDev << " "
	   "n: " << count << " "
	   "min: " << h.minimum() << " max: " << h.maximum() << '\n';
}

double PDF::getP(size_t idx) const
{
	return (idx <= m_maxIdx) ? m_dist[idx] : m_minp;
}

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
