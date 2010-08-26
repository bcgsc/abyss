#ifndef ADJIO_H
#define ADJIO_H 1

#include "ContigID.h"
#include "ContigGraph.h"
#include <cassert>
#include <istream>
#include <limits> // for numeric_limits
#include <ostream>

/** Output a contig adjacency graph. */
template <typename Graph>
std::ostream& write_adj(std::ostream& out, const Graph& g)
{
	typedef typename Graph::vertex_iterator vertex_iterator;
	typedef typename Graph::adjacency_iterator adjacency_iterator;

	std::pair<vertex_iterator, vertex_iterator> vit = g.vertices();
	bool sense = false;
	for (vertex_iterator u = vit.first; u != vit.second; ++u,
			sense = !sense) {
		if (g.is_removed(*u))
			continue;
		if (!sense)
			out << ContigID(*u) << g[*u];
		out << "\t;";
		std::pair<adjacency_iterator, adjacency_iterator>
			adj = g.adjacent_vertices(*u);
		for (adjacency_iterator v = adj.first; v != adj.second; ++v)
			out << ' ' << (*v ^ sense);
		if (sense)
			out << '\n';
	}
	return out;
}

/** Read a contig adjacency graph. */
template <typename Graph>
std::istream& read_adj(std::istream& in, ContigGraph<Graph>& g)
{
	typedef typename Graph::vertex_descriptor vertex_descriptor;
	typedef typename Graph::vertex_property_type vertex_property_type;

	// Read the vertex properties.
	g.clear();
	assert(in);
	ContigID id(-1);
	vertex_property_type prop;
	while (in >> id >> prop) {
		in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		vertex_descriptor v = g.add_vertex(prop);
		assert(v == vertex_descriptor(id, false));
	}
	assert(in.eof());
	assert(g.num_vertices() > 0);
	ContigID::lock();

	// Read the edges.
	in.clear();
	in.seekg(std::ios_base::beg);
	assert(in);
	while (in >> id) {
		in.ignore(std::numeric_limits<std::streamsize>::max(), ';');
		vertex_descriptor u(id, false);
		for (int sense = false; sense <= true; ++sense) {
			std::string s;
			getline(in, s, !sense ? ';' : '\n');
			assert(in.good());
			std::istringstream ss(s);
			for (vertex_descriptor v; ss >> v;)
				g.Graph::add_edge(u ^ sense, v ^ sense);
			assert(ss.eof());
		}
	}
	assert(in.eof());
	return in;
}

#endif
