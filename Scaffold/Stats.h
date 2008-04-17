#ifndef STATS_H
#define STATS_H

// Simple class for holding stats related to PET data
#include <vector>
#include <map>

typedef std::vector<double> PDF;
typedef std::map<int, int> histogram;

class Stats
{
	public:
		Stats();
		
		const PDF& GetPDF() const { return m_pdf; } 
		
		// Maximum Likelihood Estimator functions
		int MaxLikelihoodEst(std::vector<int>& pairDistance);
		
		// Compute the likelihood of the distribution
		double ComputeLikelihood(int d, std::vector<int>& testDist);
					
		// Generate the empirical distribution of pair distances
		void GenerateStatsFromHistogram(const histogram& h);
		
		// Generate the standard deviation of the estimate
		double GetStdDevOfEstimate(int n);	
		
		// Estimate the number of pairs that should exist between the contigs given the input parameters
		int EstimateNumberOfPairs(int estDist, int contigLength, int numberOfReads);	
				
	private:
		PDF m_pdf;
		double m_stdDev;
		double m_mean;
		int m_maxDist;
	
	
};

#endif
