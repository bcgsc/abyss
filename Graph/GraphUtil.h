#ifndef GRAPHUTIL_H
#define GRAPHUTIL_H 1

#include "Histogram.h"
#include <boost/graph/graph_traits.hpp>
#include <iomanip>
#include <ostream>

using boost::graph_traits;

/** Return the number of vertices that have been marked as removed. */
template <typename Graph>
typename graph_traits<Graph>::vertices_size_type
num_vertices_removed(const Graph& g)
{
	typedef graph_traits<Graph> GTraits;
	typedef typename GTraits::vertices_size_type vertices_size_type;
	typedef typename GTraits::vertex_iterator vertex_iterator;
	vertices_size_type n = 0;
	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u)
		if (get(vertex_removed, g, *u))
			n++;
	return n;
}

/** Print a histogram of the degree. */
template <typename Graph>
Histogram printHistogram(const Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
	Histogram h;
	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u) {
		if (get(vertex_removed, g, *u))
			continue;
		h.insert(out_degree(*u, g));
	}
	return h;
}

/** Print statistics of the number of vertices and edges. */
template <typename Graph>
std::ostream& printGraphStats(std::ostream& out, const Graph& g)
{
	using std::setprecision;
	unsigned v = num_vertices(g) - num_vertices_removed(g);
	unsigned e = num_edges(g);
	out << "V=" << v << " E=" << e
		<< " E/V=" << setprecision(3) << (float)e / v << std::endl;

	Histogram h = printHistogram(g);
	unsigned n = h.size();
	unsigned n0 = h.count(0), n1 = h.count(1), n234 = h.count(2, 4);
	unsigned n5 = n - (n0 + n1 + n234);
	return out <<
		"Degree: " << h.barplot(h.maximum() + 1) << "\n"
		"        01234\n"
		"0: " << setprecision(2) << (float)100 * n0 / n << "% "
		"1: " << setprecision(2) << (float)100 * n1 / n << "% "
		"2-4: " << setprecision(2) << (float)100 * n234 / n << "% "
		"5+: " << setprecision(2) << (float)100 * n5 / n << "% "
		"max: " << h.maximum() << std::endl;
}

/** Pass graph statistics  -- values only . */
template <typename Graph>
std::vector<int> passGraphStatsVal(const Graph& g)
{
#if _SQL
	Histogram h = printHistogram(g);
	unsigned n = h.size(),
			 n0 = h.count(0),
			 n1 = h.count(1),
			 n234 = h.count(2, 4),
			 n5 = n - (n0 + n1 + n234);

	return make_vector<int>()
		<< num_vertices(g) - num_vertices_removed(g)
		<< num_edges(g)
		<< round(100.0 * n0 / n)
		<< round(100.0 * n1 / n)
		<< round(100.0 * n234 / n)
		<< round(100.0 * n5 / n)
		<< h.maximum();
#else
	(void)g;
	return make_vector<int>();
#endif
}
#endif
