#ifndef GRAPHUTIL_H
#define GRAPHUTIL_H 1

#include "Histogram.h"
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
	return out <<
		"Degree: " << h.barplot(h.maximum() + 1) << "\n"
		"        01234\n";
}

#endif
