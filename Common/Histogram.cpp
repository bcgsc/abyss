#include "Histogram.h"

/** Remove samples less than the specified threshold. */
Histogram Histogram::trimLow(T threshold) const
{
	Histogram h;
	for (Histogram::Map::const_iterator it = m_map.begin();
			it != m_map.end(); it++)
		if (it->first >= threshold)
			h.insert(it->first, it->second);
	return h;
}

/** Trim off the bottom fraction/2 and top fraction/2 data points.
 * At least (1 - fraction) of the data will remain.
 */
Histogram Histogram::trimFraction(double fraction) const
{
	double low_cutoff = fraction/2;
	double high_cutoff = 1.0f - fraction/2;
	unsigned n = size();

	double cumulative = 0;
	Histogram newHist;
	for (Histogram::Map::const_iterator it = m_map.begin();
			it != m_map.end(); it++) {
		double temp_total = cumulative + (double)it->second / n;
		if (temp_total > low_cutoff && cumulative < high_cutoff)
			newHist.insert(it->first, it->second);
		cumulative = temp_total;
	}

	return newHist;
}
