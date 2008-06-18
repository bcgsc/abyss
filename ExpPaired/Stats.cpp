#include "Stats.h"
#include <math.h>
#include <iostream>

const double MINP = 0.00001f;

// Functions

void KLDiv(const PDF& p, const PDF& q)
{
	p.print();
	q.print();
	
	const size_t maxIdx = q.getMax();
	
	double sum = 0.0f;
	
	printf("KL:\n");
	for(size_t idx = 0; idx <= maxIdx; ++idx)
	{
		double val = p.getP(idx) * log(p.getP(idx) / q.getP(idx));
		sum += val;
		printf("%zu p: %lf q: %lf %lf (%lf)\n", idx, p.getP(idx), q.getP(idx), val, sum);
	}
	
	printf("KLDiv: %lf\n", sum);
}

void ChiSquare(const PDF& ref, const Histogram& sample)
{
	// Count the number of items in the sample hist
	int count = sample.getSumCount();
	count /= 6;
	size_t maxIdx = ref.getMax();
	
	// Compute the expected vector
	double chi = 0.0f;
	for(size_t idx = 0; idx <= maxIdx; ++idx)
	{
		double expected = ref.getP(idx) * count;
		double observed = (sample.getCount(idx) / 6);
		double val = pow((expected - observed), 2.0f) / expected;
		printf("idx %zu Observed %lf Expect %lf val %lf\n", idx, observed, expected, val);
		chi += val;
	}
	
	printf("Chi-sq: %lf\n", chi);
	
}

// Classes
void Histogram::addDataPoint(int data)
{
	m_data[data]++;
}

int Histogram::getMin() const
{
	if(m_data.empty())
	{
		return 0;
	}
	else
	{
		return m_data.begin()->first;
	}	
}

int Histogram::getMax() const
{
	if(m_data.empty())
	{
		return 0;
	}
	else
	{
		return m_data.rbegin()->first;
	}
}

int Histogram::getCount(int index) const
{
	IntIntMap::const_iterator iter = m_data.find(index);
	if(m_data.find(index) == m_data.end())
	{
		return 0;
	}
	else
	{
		return iter->second;
	}
}

int Histogram::getSumCount() const
{
	int min = 0;
	int max = getMax();
	
	int sum = 0;
	for(int idx = min; idx <= max; ++idx)
	{
		sum += getCount(idx);
	}
	return sum;
}

void Histogram::print() const
{
	
	int min = 0;
	int max = getMax();
	
	printf("Hist: [%d-%d]\n", min, max);
	for(int idx = min; idx <= max; ++idx)
	{
		printf("%d %d\n", idx, getCount(idx));
	}
}

// Construct a pdf from a histogram
PDF::PDF(const Histogram& h)
{
	m_maxVal = h.getMax();
	double count = 0;

	for(IntIntMap::const_iterator histIter = h.m_data.begin(); histIter != h.m_data.end(); histIter++)
	{
		count += histIter->second;
	}
	
	// Create the initial pdf with all values being 0
	m_dist = DoubleVec(m_maxVal+1, 0.0f);
	for(size_t i = 0; i <= m_maxVal; i++)
	{
		int v = h.getCount(i);

		if(v > 0)
		{
			m_dist[i] = static_cast<double>(v) / count;
		}
		else
		{
			m_dist[i] = 0;
		}
	}	
}

double PDF::getP(size_t idx) const 
{ 
	return (idx <= m_maxVal) ? m_dist[idx] : MINP; 
}

//
//
//
void PDF::print() const
{
	for(size_t idx = 0; idx <= m_maxVal; ++idx)
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
	for(size_t idx = 0; idx <= m_maxVal; ++idx)
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
		size_t nextRight = (rightIdx < m_maxVal) ? rightIdx + 1 : rightIdx;
		
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

PairedStats::PairedStats() : m_maxDist(0)
{
	
}

//
// Perform a maximum likelihood estimate over the pdf and input distribution
//
int PairedStats::maxLikelihoodEst(int min, int max, std::vector<int>& pairDistance)
{
	double maxL = -999999;
	int bestDist = min;
	
	int minDist = min;
	int maxDist = max;
	for(int i = minDist; i < maxDist; i++)
	{
		double v = 	computeLikelihood(i, pairDistance);
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
double PairedStats::computeLikelihood(int d, std::vector<int>& testDist)
{
	double sum = 0.0f;

	for(std::vector<int>::iterator iter = testDist.begin(); iter != testDist.end(); iter++)
	{
		int val = *iter + d;
		double p;
		
		// Is this value in range of the pdf?
		if(val >= 0 && val < (int)m_pdf.m_dist.size())
		{
			p = m_pdf.m_dist[val];	
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
int PairedStats::estimateNumberOfPairs(size_t estDist, int contigLength, int numberOfReads)
{
	double p = 0.0f;
	for(size_t i = estDist; i < estDist + contigLength && i < m_maxDist; i++)
	{
		p += m_pdf.m_dist[i];
	}
	
	assert(p <= 1.0f);
	return (int)round(p * (double)numberOfReads);
}

//
// Convert the histogram to a probability density function
void PairedStats::generateStats(const Histogram& h)
{
	m_maxDist = h.getMax();
	double count = 0;
	double sum = 0;
	for(IntIntMap::const_iterator histIter = h.m_data.begin(); histIter != h.m_data.end(); histIter++)
	{
		count += histIter->second;
		sum += (histIter->second * histIter->first);
	}
	
	double mean = sum / count;
	
	// Compute the standard devition
	sum = 0;
	count = 0;
	for(IntIntMap::const_iterator histIter = h.m_data.begin(); histIter != h.m_data.end(); histIter++)
	{

		count += histIter->second;
		sum += (static_cast<double>(histIter->second) * pow((static_cast<double>(histIter->first) - mean), 2));
	}
	
	sum /= count;
	m_stdDev = sqrt(sum);
	m_mean = mean;
	printf("mean distance: %lf std dev: %lf count: %d\n", m_mean, m_stdDev, (int)count);
	
	// Create the initial pdf with all values being 0
	m_pdf = PDF(h);
}

//
//
//
void PairedStats::computeCDF(int low, int high)
{
	double sum = 0.0f;
	
	size_t start = (low > 0) ? low : 0;
	size_t end = (high < (int)m_maxDist) ? high : m_maxDist;
	
	for(size_t idx = start; idx < end; ++idx)
	{	
		sum += m_pdf.m_dist[idx];
	}
	
	printf("CumulativeProb [%zu-%zu] = %lf\n", start, end, sum);
}

//
// Get the standard deviation of the estimate based on the empirical standard deviation and the number of data points
//
double PairedStats::getStdDevOfEstimate(int n)
{
	return m_stdDev / sqrt((double)n);
}
