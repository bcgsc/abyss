#include "MLE.h"
#include "PDF.h"
#include <algorithm> // for swap
#include <cassert>
#include <climits> // for INT_MAX
#include <limits> // for numeric_limits

using namespace std;

namespace opt {
	extern unsigned k;
}

/** This window function is a triangle with a flat top, or a rectangle
 * with sloped sides.
 */
class WindowFunction {
	public:
		WindowFunction(int len0, int len1, int d)
		{
			assert(len0 > 0);
			assert(len1 > 0);
			assert(len0 <= len1);
			if (d < 0)
				d = 0;
			x0 = d;
			x1 = d + len0;
			x2 = d + len1;
			x3 = d + len0 + len1;
			y = len0;
		}

		/** Return this window function evaluated at x. */
		double operator ()(int x) const
		{
			return (x <= x0 ? 1
					: x < x1 ? x - x0
					: x < x2 ? y
					: x < x3 ? x3 - x
					: 1) / (double)y;
		}

	private:
		/** Parameters of this window function. */
		int x0, x1, x2, x3, y;
};

/** Compute the log likelihood that these samples came from the
 * specified distribution shifted by the parameter theta and scaled by
 * the specified window function.
 * @param theta the parameter of the PDF
 * @param samples the samples
 * @param pdf the PDF scaled by the window function
 * @param window the window function used to scale the PDF
 * @param n [out] the number of samples with a non-zero probability
 * @return the log likelihood
 */
static double computeLikelihood(int theta, const vector<int>& samples,
		const PDF& pdf, const WindowFunction& window,
		unsigned &n)
{
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
		sum += log(p * window(val));
		if (p > pdf.getMinP())
			n++;
	}
	return sum;
}

/** Return the most likely distance between two contigs. */
static int maximumLikelihoodEstimate(int first, int last,
		const vector<int>& samples,
		const PDF& pdf, const WindowFunction& window,
		unsigned &n)
{
	double bestLikelihood = -numeric_limits<double>::max();
	int bestTheta = first;
	unsigned bestn = 0;
	for (int theta = first; theta < last; theta++) {
		unsigned trialn;
		double likelihood = computeLikelihood(theta, samples,
				pdf, window, trialn);
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
 * @param len0 the length of the first contig in bp
 * @param len1 the length of the second contig in bp
 */
int maximumLikelihoodEstimate(int first, int last,
		const vector<int>& samples, const PDF& pdf,
		unsigned len0, unsigned len1,
		unsigned& n)
{
	// Convert the lengths of the contigs from bp to k-mer.
	len0 -= opt::k - 1;
	len1 -= opt::k - 1;
	if (len0 > len1)
		swap(len0, len1);
	int d0 = maximumLikelihoodEstimate(first, last, samples,
			pdf, WindowFunction(len0, INT_MAX/2, 0), n);
	return maximumLikelihoodEstimate(first, last, samples,
			pdf, WindowFunction(len0, len1, d0), n);
}
