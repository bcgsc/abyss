#include "MLE.h"
#include "Stats.h"
#include <cassert>
#include <climits> // for INT_MAX
#include <limits> // for numeric_limits

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
	double bestLikelihood = -numeric_limits<double>::max();
	int bestTheta = first;
	unsigned bestn = 0;
	for (int theta = first; theta < last; theta++) {
		unsigned trialn;
		double likelihood = computeLikelihood(theta, samples, pdf,
				len0, len1, d, trialn);
		if (likelihood > bestLikelihood) {
			bestLikelihood = likelihood;
			bestTheta = theta;
			bestn = trialn;
		}
	}
	n = bestn;
	return bestTheta;
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
				INT_MAX/4, INT_MAX/4, 0, n);
}
