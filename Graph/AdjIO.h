#ifndef ADJIO_H
#define ADJIO_H 1

#include "ContigID.h"
#include "ContigGraph.h"
#include <cassert>
#include <istream>
#include <limits> // for numeric_limits
#include <ostream>

template <typename EdgeProp>
void write_edge_prop(std::ostream& out, const EdgeProp& ep)
{
	if (ep != EdgeProp())
		out << " [" << ep << ']';
}

static inline void write_edge_prop(std::ostream&, const no_property&)
{
}

/** Output a contig adjacency graph. */
template <typename Graph>
std::ostream& write_adj(std::ostream& out, const Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	typedef typename graph_traits<Graph>::vertex_iterator
		vertex_iterator;
	typedef typename graph_traits<Graph>::out_edge_iterator
		out_edge_iterator;

	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	bool sense = false;
	for (vertex_iterator u = vit.first; u != vit.second; ++u,
			sense = !sense) {
 		if (get(vertex_removed, g, *u))
			continue;
		if (!sense)
			out << ContigID(*u) << get(vertex_bundle, g, *u);
		out << "\t;";
		std::pair<out_edge_iterator, out_edge_iterator>
			adj = out_edges(*u, g);
		for (out_edge_iterator e = adj.first; e != adj.second; ++e) {
			vertex_descriptor v = target(*e, g);
			assert(!get(vertex_removed, g, v));
			out << ' ' << (v ^ sense);
			write_edge_prop(out, get(edge_bundle, g, e));
		}
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
	typedef typename Graph::edge_property_type edge_property_type;

	// Read the vertex properties.
	g.clear();
	assert(in);
	ContigID id(-1);
	vertex_property_type prop;
	while (in >> id >> prop) {
		in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		vertex_descriptor v = add_vertex(prop, g);
		assert(v == vertex_descriptor(id, false));
	}
	assert(in.eof());
	assert(num_vertices(g) > 0);
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
			for (vertex_descriptor v; ss >> v >> std::ws;) {
				if (ss.peek() == '[') {
					ss.get();
					edge_property_type ep;
					ss >> ep;
					ss.ignore(std::numeric_limits<
							std::streamsize>::max(), ']');
					g.Graph::add_edge(u ^ sense, v ^ sense, ep);
				} else
					g.Graph::add_edge(u ^ sense, v ^ sense);
			}
			assert(ss.eof());
		}
	}
	assert(in.eof());
	return in;
}

template <typename Graph>
class AdjWriter
{
	const Graph& m_g;
  public:
	AdjWriter(const Graph& g) : m_g(g) { }
	friend std::ostream& operator<<(std::ostream& out,
			const AdjWriter& o)
	{
		return write_adj<Graph>(out, o.m_g);
	}
};

template <typename Graph>
AdjWriter<Graph> adj_writer(const Graph& g)
{
	return AdjWriter<Graph>(g);
}

#endif
