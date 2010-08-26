#ifndef DOTIO_H
#define DOTIO_H 1

#include "Graph.h"
#include <ostream>

/** Output a GraphViz dot graph. */
template <typename Graph>
std::ostream& write_dot(std::ostream& out, const Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_iterator
		vertex_iterator;
	typedef typename graph_traits<Graph>::adjacency_iterator
		adjacency_iterator;
	typedef typename vertex_property<Graph>::type
		vertex_property_type;

	std::pair<vertex_iterator, vertex_iterator>
		vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u) {
		if (is_removed(*u, g))
			continue;
		out << vertex_property_writer<Graph, vertex_property_type>(
					g, *u);
		unsigned outdeg = out_degree(*u, g);
		if (outdeg == 0)
			continue;
		out << '"' << *u << "\" ->";
		if (outdeg > 1)
			out << " {";
		std::pair<adjacency_iterator, adjacency_iterator>
			adj = adjacent_vertices(*u, g);
		for (adjacency_iterator v = adj.first; v != adj.second; ++v)
			out << " \"" << *v << '"';
		if (outdeg > 1)
			out << " }";
		out << '\n';
	}
	return out;
}

/** Output a GraphViz dot graph. */
template <typename Graph>
struct dot_writer
{
	const Graph& g;
	dot_writer(const Graph& g) : g(g) { }
	friend std::ostream& operator<<(std::ostream& out,
			const dot_writer& o)
	{
		return write_dot<Graph>(out, o.g);
	}
};

#endif
