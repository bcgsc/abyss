#include "Histogram.h"
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <string>

using namespace std;

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
	size_type n = size();

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

/** Bin these elements into n buckets. The number of buckets returned
 * may be smaller than n.
 */
Histogram::Bins Histogram::bin(unsigned n) const
{
	Histogram::Bins bins;
	bins.reserve(n);
	T nperbucket = (T)ceilf((float)(maximum() - minimum()) / n);
	T next = minimum() + nperbucket;
	Histogram::Bins::value_type count = 0;
	for (Histogram::Map::const_iterator it = m_map.begin();
			it != m_map.end(); it++) {
		if (it->first >= next) {
			bins.push_back(count);
			count = 0;
			next += nperbucket;
		}
		count += it->second;
	}
	if (count > 0)
		bins.push_back(count);
	return bins;
}

/** Return a unicode bar plot.
 * @param n number of buckets
 */
string Histogram::barplot(unsigned nbins) const
{
	/** Unicode bar characters. */
	static const char* bars[10] = {
		" ", "_",
		"\342\226\201", // 9601
		"\342\226\202", // 9602
		"\342\226\203", // 9603
		"\342\226\204", // 9604
		"\342\226\205", // 9605
		"\342\226\206", // 9606
		"\342\226\207", // 9607
		"\342\226\210", // 9608
	};

	Histogram::Bins bins = bin(nbins);
	ostringstream ss;
	Histogram::Bins::value_type max
		= 1 + *max_element(bins.begin(), bins.end());
	for (Histogram::Bins::const_iterator it = bins.begin();
			it != bins.end(); ++it)
		ss << bars[10 * *it / max];
	string s(ss.str());
	while (!s.empty() && *s.rbegin() == ' ')
		s.erase(s.size() - 1);
	return s;
}

/** Return a unicode bar plot. */
string Histogram::barplot() const
{
	const char *columns = getenv("COLUMNS");
	return barplot(columns == NULL ? 80 : strtoul(columns, NULL, 0));
}
