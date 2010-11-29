#ifndef GRAPHUTIL_H
#define GRAPHUTIL_H 1

#include "Histogram.h"
#include <iomanip>
#include <ostream>

template <typename Graph>
std::ostream& printGraphStats(std::ostream& out, const Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_iterator
		vertex_iterator;

	unsigned v = num_vertices(g);
	unsigned e = num_edges(g);
	out << "V=" << v << " E=" << e
		<< " E/V=" << (float)e / v << std::endl;

	// Print a histogram of the degree.
	Histogram h;
	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u)
		h.insert(out_degree(*u, g));
	unsigned n = h.size();
	unsigned n0 = h.count(0), n1 = h.count(1), n234 = h.count(2, 4);
	unsigned n5 = n - (n0 + n1 + n234);
	using std::setprecision;
	return out <<
		"Degree: " << h.barplot(h.maximum() + 1) << "\n"
		"        01234\n"
		"0: " << setprecision(2) << (float)100 * n0 / n << "% "
		"1: " << setprecision(2) << (float)100 * n1 / n << "% "
		"2-4: " << setprecision(2) << (float)100 * n234 / n << "% "
		"5+: " << setprecision(2) << (float)100 * n5 / n << "% "
		"max: " << h.maximum() << std::endl;
}

#endif
