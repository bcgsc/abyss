#ifndef DOTIO_H
#define DOTIO_H 1

#include "Graph.h"
#include <ostream>

template <typename V, typename VertexProp>
std::ostream& write_vertex(std::ostream& out,
		V u, const VertexProp& vp)
{
	return out << '"' << u << "\" [" << vp << "]\n";
}

template <typename V>
std::ostream& write_vertex(std::ostream& out, V, no_property)
{
	return out;
}

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
 		if (get(vertex_removed, g, *u))
			continue;
		write_vertex(out, *u, get(vertex_bundle, g, *u));
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
struct DotWriter
{
	const Graph& g;
	DotWriter(const Graph& g) : g(g) { }
	friend std::ostream& operator<<(std::ostream& out,
			const DotWriter& o)
	{
		return write_dot<Graph>(out, o.g);
	}
};

/** Output a GraphViz dot graph. */
template <typename Graph>
DotWriter<Graph> dot_writer(const Graph& g)
{
	return DotWriter<Graph>(g);
}

#endif
