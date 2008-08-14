#ifndef STATS_H
#define STATS_H

// Simple class for holding stats related to PET data
#include <vector>
#include <map>
#include <math.h>

// Classes
typedef std::map<int, int> IntIntMap;

struct Histogram
{
	Histogram() {}
	Histogram(const std::vector<int>& data);
	
	void addDataPoint(int data);
	void addMultiplePoints(int value, int count);
	
	Histogram trim(double percent);
	
	int getSumCount() const;
	int getCount(int index) const;
	int getMin() const;
	int getMax() const;
	
	void print() const;
	
	IntIntMap m_data;
};

typedef std::vector<double> DoubleVec;
struct PDF
{
	PDF() {};
	PDF(const Histogram& h);
	
	double getP(size_t idx) const;
	
	size_t getMaxIdx() const { return m_maxIdx; }
	void print() const;
	int getSampleStdDev(int n) const { return (int)ceil(m_stdDev / sqrt((double)n)); }
	
	size_t m_maxIdx;
	DoubleVec m_dist;
	double m_mean;
	double m_stdDev;
	
	// calculate the minimal range in which p% of the values will fall into
	void calculateMinimalRange(double p, size_t& low, size_t& high) const;
};

struct CDF
{
	CDF() {};
	CDF(const PDF& pdf);
	
	double getP(size_t idx) const;
	
	size_t getMaxIdx() const { return m_maxIdx; }
	void print() const;
	
	size_t m_maxIdx;
	DoubleVec m_dist;
};


// Functions
void KLDiv(const PDF& p, const PDF& q);
void ChiSquare(const PDF& ref, const Histogram& sample);
bool KSTestCont(std::vector<int> observations, const PDF& p);
double approximateKSCritValue(int n, double alpha);

// Maximum Likelihood Estimator functions
int maxLikelihoodEst(int min, int max, std::vector<int>& pairDistance, const PDF& pdf, double& ratio);

// Compute the likelihood of the distribution
double computeLikelihood(int param, std::vector<int>& testDist, const PDF& pdf);

class PairedStats
{
	public:
		PairedStats();
		
		const PDF& getPDF() const { return m_pdf; } 
					
		// Generate stats from the histogram
		void generateStats(const Histogram& h);
		
		// Generate the standard deviation of the estimate
		double getStdDevOfEstimate(int n);	
		
		// Estimate the number of pairs that should exist between the contigs given the input parameters
		int estimateNumberOfPairs(size_t estDist, int contigLength, int numberOfReads);

		// Compute the sum of probablities between low and high
		void computeCDF(int low, int high);
		
		// Return the maximum distance
		size_t getMax() const { return m_maxDist; }
		
		double getMean() const { return m_mean; }
		double getStdDev() const { return m_stdDev; }
				
	private:
		PDF m_pdf;
		double m_stdDev;
		double m_mean;
		size_t m_maxDist;
	
	
};

#endif
