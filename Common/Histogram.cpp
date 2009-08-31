#include "Histogram.h"

/** Trim off the bottom percent/2 and top percent/2 data points.
 * At least (1 - percent) of the data will remain.
 */
Histogram Histogram::trim(double percent) const
{
	double low_cutoff = percent/2;
	double high_cutoff = 1.0f - percent/2;
	unsigned n = size();

	double cumulative = 0;
	Histogram newHist;
	for (Histogram::Map::const_iterator it = begin();
			it != end(); it++) {
		double temp_total = cumulative + (double)it->second / n;
		if (temp_total > low_cutoff && cumulative < high_cutoff)
			newHist.insert(it->first, it->second);
		cumulative = temp_total;
	}

	return newHist;
}
