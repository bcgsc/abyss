#ifndef MLE_H
#define MLE_H 1

#include <vector>

class PDF;

int maximumLikelihoodEstimate(int first, int last,
		const std::vector<int>& samples, const PDF& pdf,
		unsigned len0, unsigned len1, bool rf, unsigned& n);

#endif
