#ifndef STATS_H
#define STATS_H

// Simple class for holding stats related to PET data
#include <vector>
#include <map>

// Classes
typedef std::map<int, int> IntIntMap;

struct Histogram
{
	Histogram() {}
	
	void addDataPoint(int data);
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
	
	size_t getMax() const { return m_maxVal; }
	void print() const;
	
	size_t m_maxVal;
	DoubleVec m_dist;
	
	// calculate the minimal range in which p% of the values will fall into
	void calculateMinimalRange(double p, size_t& low, size_t& high) const;
};


// Functions
void KLDiv(const PDF& p, const PDF& q);
void ChiSquare(const PDF& ref, const Histogram& sample);

class PairedStats
{
	public:
		PairedStats();
		
		const PDF& getPDF() const { return m_pdf; } 
		
		// Maximum Likelihood Estimator functions
		int maxLikelihoodEst(int min, int max, std::vector<int>& pairDistance);
		
		// Compute the likelihood of the distribution
		double computeLikelihood(int d, std::vector<int>& testDist);
					
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
