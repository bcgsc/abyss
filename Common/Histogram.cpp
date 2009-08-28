#include "Histogram.h"

/** Trim off the bottom percent/2 and top percent/2 data points.
 * At least (1 - percent) of the data will remain.
 */
Histogram Histogram::trim(double percent)
{
	double half_percent = percent/2;
	double low_cutoff = half_percent;
	double high_cutoff = 1.0f - half_percent;
	double total = (double)getSumCount();

	T min = 0;
	T max = getMax();
	double cumulative = 0.0f;

	Histogram newHist;
	for (T value = min; value <= max; ++value) {
		unsigned currCount = getCount(value);
		double frac = currCount / total;
		double temp_total = cumulative + frac;
		if (temp_total > low_cutoff && cumulative < high_cutoff)
			newHist.addMultiplePoints(value, currCount);
		cumulative = temp_total;
	}

	return newHist;
}
