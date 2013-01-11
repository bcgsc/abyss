#include "MLE.h"
#include "PMF.h"
#include <boost/tuple/tuple.hpp>
#include <algorithm> // for swap
#include <cassert>
#include <limits> // for numeric_limits
#include <utility>

using namespace std;
using boost::tie;

/** This window function is a triangle with a flat top, or a rectangle
 * with sloped sides.
 */
class WindowFunction {
	public:
		WindowFunction(int len0, int len1)
			: x1(len0), x2(len1), x3(len0 + len1)
		{
			assert(len0 > 0);
			assert(len1 > 0);
			assert(len0 <= len1);
		}

		/** Return this window function evaluated at x. */
		double operator()(int x) const
		{
			return (x <= 0 ? 1
					: x < x1 ? x
					: x < x2 ? x1
					: x < x3 ? x3 - x
					: 1) / (double)x1;
		}

	private:
		/** Parameters of this window function. */
		int x1, x2, x3;
};

/** Compute the log likelihood that these samples came from the
 * specified distribution shifted by the parameter theta.
 * @param theta the parameter of the PMF, f_theta(x)
 * @param samples the samples
 * @param pmf the probability mass function
 * @return the log likelihood
 */
static pair<double, unsigned>
computeLikelihood(int theta, const Histogram& samples, const PMF& pmf)
{
	double likelihood = 0;
	unsigned nsamples = 0;
	for (Histogram::const_iterator it = samples.begin();
			it != samples.end(); ++it) {
		double p = pmf[it->first + theta];
		unsigned n = it->second;
		likelihood += n * log(p);
		if (p > pmf.minProbability())
			nsamples += n;
	}
	return make_pair(likelihood, nsamples);
}

/** Return the most likely distance between two contigs and the number
 * of pairs that support that estimate. */
static pair<int, unsigned>
maximumLikelihoodEstimate(int first, int last,
		const Histogram& samples,
		const PMF& pmf,
		unsigned len0, unsigned len1)
{
	first = max(first, (int)pmf.minValue() - samples.maximum());
	last = min(last, (int)pmf.maxValue() - samples.minimum());

	/* When randomly selecting fragments that span a given point,
	 * longer fragments are more likely to be selected than
	 * shorter fragments.
	 */
	WindowFunction window(len0, len1);

	unsigned nsamples = samples.size();
	double bestLikelihood = -numeric_limits<double>::max();
	int bestTheta = first;
	unsigned bestn = 0;
	for (int theta = first; theta <= last; theta++) {
		// Calculate the normalizing constant of the PMF, f_theta(x).
		double c = 0;
		for (int i = pmf.minValue(); i <= (int)pmf.maxValue(); ++i)
			c += pmf[i] * window(i - theta);

		double likelihood;
		unsigned n;
	   	tie(likelihood, n) = computeLikelihood(theta, samples, pmf);
		likelihood -= nsamples * log(c);
		if (n > 0 && likelihood > bestLikelihood) {
			bestLikelihood = likelihood;
			bestTheta = theta;
			bestn = n;
		}
	}
	return make_pair(bestTheta, bestn);
}

/** Return the most likely distance between two contigs and the number
 * of pairs that support that distance estimate.
 * @param len0 the length of the first contig in bp
 * @param len1 the length of the second contig in bp
 * @param rf whether the fragment library is oriented reverse-forward
 * @param[out] n the number of samples with a non-zero probability
 */
int maximumLikelihoodEstimate(unsigned l,
		int first, int last,
		const vector<int>& samples, const PMF& pmf,
		unsigned len0, unsigned len1, bool rf,
		unsigned& n)
{
	assert(first < last);
	assert(!samples.empty());

	// The aligner is unable to map reads to the ends of the sequence.
	// Correct for this lack of sensitivity by subtracting l-1 bp from
	// the length of each sequence, where the aligner requires a match
	// of at least l bp. When the fragment library is oriented
	// forward-reverse, subtract 2*(l-1) from each sample.
	assert(l > 0);
	assert(len0 >= l);
	assert(len1 >= l);
	len0 -= l - 1;
	len1 -= l - 1;

	if (len0 > len1)
		swap(len0, len1);

	if (rf) {
		// This library is oriented reverse-forward.
		Histogram h(samples.begin(), samples.end());
		int d;
		tie(d, n) = maximumLikelihoodEstimate(
				first, last, h,
				pmf, len0, len1);
		return d;
	} else {
		// This library is oriented forward-reverse.
		// Subtract 2*(l-1) from each sample.
		Histogram h;
		typedef vector<int> Samples;
		for (Samples::const_iterator it = samples.begin();
				it != samples.end(); ++it) {
			assert(*it > 2 * (int)(l - 1));
			h.insert(*it - 2 * (l - 1));
		}
		int d;
		tie(d, n) = maximumLikelihoodEstimate(
				first, last, h,
				pmf, len0, len1);
		return max(first, d - 2 * (int)(l - 1));
	}
}
