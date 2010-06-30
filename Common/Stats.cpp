#include "Stats.h"
#include "Histogram.h"
#include <iomanip>
#include <iostream>

using namespace std;

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
