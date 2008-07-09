#include "Stats.h"
#include <math.h>
#include <iostream>
#include <algorithm>

const double MINP = 0.00001f;

// Functions

void KLDiv(const PDF& p, const PDF& q)
{
	p.print();
	q.print();
	
	const size_t maxIdx = q.getMaxIdx();
	
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

void KSTestCont(std::vector<int> observations, const PDF& p)
{
	// convert the observations vector to a histogram
	Histogram hist(observations);

	// convert the histogram to a pdf (a frequency distribution in this case)
	PDF freqs(hist);
	
	// Convert both distributions to CDFs
	CDF obsCDF(freqs);
	CDF refCDF(p);
	
	// calculate the test statistic which is max[ max(Di), max(Di') ]
	// where Di = abs(Oi - Fi)
	// and Di' = abs(Oi-1 - Fi)
	size_t maxIdx = std::max(refCDF.getMaxIdx(), obsCDF.getMaxIdx());
	double D = 0.f;
	for(size_t idx = 0; idx <= maxIdx; ++idx)
	{
		double Di = fabs(obsCDF.getP(idx) - refCDF.getP(idx));
		double Diprev;
		if(idx == 0)
		{
			Diprev = refCDF.getP(idx);
		}
		else
		{
			Diprev = fabs(obsCDF.getP(idx - 1) - refCDF.getP(idx));
		}
		
		D = std::max(D, Di);
		D = std::max(D, Diprev);
		
		printf("cdf %lf %lf %lf %lf %lf\n", obsCDF.getP(idx), refCDF.getP(idx), Di, Diprev, D);
	}
	
	double alpha = 0.01;
	double crit = approximateKSCritValue(observations.size(), alpha);
	printf("D: %lf Crit: %lf Result: %d\n", D, crit, D < crit);
}

double approximateKSCritValue(int n, double alpha)
{
	// Formula from JH Zar - Biostatistical Analysis, 2nd ed, pg 548
	// The approximation is accurate to ~3.5% for n=10, alpha = 0.001
	double crit = sqrt(-log(alpha/2) / (2*n)) - (0.16693/n);
	return crit;
}

void ChiSquare(const PDF& ref, const Histogram& sample)
{
	// Count the number of items in the sample hist
	int count = sample.getSumCount();
	count /= 6;
	size_t maxIdx = ref.getMaxIdx();
	
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

//
// Perform a maximum likelihood estimate over the pdf and input distribution
//
int maxLikelihoodEst(int min, int max, std::vector<int>& pairDistance, const PDF& pdf, double& ratio)
{
	double maxL = -999999;
	double nextBestL = maxL;
	int bestDist = min;
	
	int minDist = min;
	int maxDist = max;
	for(int i = minDist; i < maxDist; i++)
	{
		double v = computeLikelihood(i, pairDistance, pdf);
		if(v > maxL)
		{
			nextBestL = maxL;
			maxL = v;
			bestDist = i;
		}
	}
	
	// Compute the ratio between the best and second best
	ratio = maxL / nextBestL;
	
	return bestDist;
}

//
// Compute the log likelihood function over the test distribution
//
double computeLikelihood(int param, std::vector<int>& testDist, const PDF& pdf)
{
	double sum = 0.0f;

	for(std::vector<int>::iterator iter = testDist.begin(); iter != testDist.end(); iter++)
	{
		int val = *iter + param;
		double p = pdf.getP(val);
		sum += log(p);
	}
	
	return sum;
}

// Classes
Histogram::Histogram(const std::vector<int>& data)
{
	for(std::vector<int>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
	{
		addDataPoint(*iter);
	}
}

void Histogram::addDataPoint(int data)
{
	m_data[data]++;
}

void Histogram::addMultiplePoints(int value, int count)
{
	m_data[value] += count;
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
	m_maxIdx = h.getMax();
	double count = 0;

	for(IntIntMap::const_iterator histIter = h.m_data.begin(); histIter != h.m_data.end(); histIter++)
	{
		count += histIter->second;
	}
	
	// Create the initial pdf with all values being 0
	m_dist = DoubleVec(m_maxIdx+1, 0.0f);
	for(size_t i = 0; i <= m_maxIdx; i++)
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
	
	// Calculate the mean
	double sum = 0;
	double sumsqr = 0;
	for(IntIntMap::const_iterator histIter = h.m_data.begin(); histIter != h.m_data.end(); histIter++)
	{
		sum += (histIter->second * histIter->first);
		sumsqr += histIter->second * (pow(histIter->first, 2.0f));
	}
	
	m_mean = sum / count;
	
	double msqr = pow(sum,2.0) / count;
	double ss1 = sumsqr - msqr;
	double t1 = ss1 / count;

	m_stdDev = sqrt(t1);
	
	printf("Calculated stats - mean: %lf stddev: %lf count: %lf\n", m_mean, m_stdDev, count);
}

double PDF::getP(size_t idx) const 
{ 
	return (idx <= m_maxIdx) ? m_dist[idx] : MINP; 
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

// Construct a cdf from a pdf
CDF::CDF(const PDF& pdf)
{
	m_maxIdx = pdf.getMaxIdx();
	double sum = 0.0f;
	// Create the initial pdf with all values being 0
	m_dist = DoubleVec(m_maxIdx+1, 0.0f);
	for(size_t i = 0; i <= m_maxIdx; i++)
	{
		sum += pdf.getP(i);
		m_dist[i] = sum;
	}	
}

double CDF::getP(size_t idx) const 
{ 
	return (idx <= m_maxIdx) ? m_dist[idx] : m_dist[m_maxIdx]; 
}

//
//
//
void CDF::print() const
{
	for(size_t idx = 0; idx <= m_maxIdx; ++idx)
	{
		printf("%zu %lf\n", idx, m_dist[idx]);
	}
}

PairedStats::PairedStats() : m_maxDist(0)
{
	
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
	assert(false && "fix std dev calc");
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
