#ifndef MLE_H
#define MLE_H 1

#include <vector>

class PMF;

int maximumLikelihoodEstimate(int first, int last,
		const std::vector<int>& samples, const PMF& pmf,
		unsigned len0, unsigned len1, bool rf, unsigned& n);

#endif
