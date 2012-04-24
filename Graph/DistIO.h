#ifndef DISTIO_H
#define DISTIO_H 1

#include <boost/graph/graph_traits.hpp>
#include <cassert>
#include <ostream>

using boost::graph_traits;

/** Output a distance estimate graph. */
template <typename Graph>
std::ostream& write_dist(std::ostream& out, const Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::vertex_iterator Vit;
	typedef typename graph_traits<Graph>::out_edge_iterator Eit;

	std::pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
 		if (get(vertex_removed, g, u)
				|| out_degree(u, g) + in_degree(u, g) == 0)
			continue;
		bool sense = get(vertex_sense, g, u);
		if (!sense)
			out << get(vertex_contig_name, g, u);
		else
			out << " ;";
		std::pair<Eit, Eit> erange = out_edges(u, g);
		for (Eit eit = erange.first; eit != erange.second; ++eit) {
			V v = target(*eit, g) ^ sense;
			assert(!get(vertex_removed, g, v));
			out << ' ' << get(vertex_name, g, v)
				<< ',' << get(edge_bundle, g, eit);
		}
		if (sense)
			out << '\n';
	}
	return out;
}

#endif
