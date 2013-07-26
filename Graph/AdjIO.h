#ifndef ADJIO_H
#define ADJIO_H 1

#include "ContigID.h"
#include "ContigGraph.h"
#include "IOUtil.h"
#include <boost/graph/graph_traits.hpp>
#include <algorithm> // for count
#include <cassert>
#include <cstdlib>
#include <istream>
#include <ostream>
#include <utility>

using namespace std::rel_ops;
using boost::graph_traits;

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
			out << get(vertex_contig_name, g, *u)
				<< get(vertex_bundle, g, *u);
		out << "\t;";
		std::pair<out_edge_iterator, out_edge_iterator>
			adj = out_edges(*u, g);
		for (out_edge_iterator e = adj.first; e != adj.second; ++e) {
			vertex_descriptor v = target(*e, g);
			assert(!get(vertex_removed, g, v));
			out << ' ' << get(vertex_name, g, v ^ sense);
			write_edge_prop(out, get(edge_bundle, g, e));
		}
		if (sense)
			out << '\n';
	}
	return out;
}

/** Read the edges of a graph in dist format. */
template <typename Graph, typename BetterEP>
std::istream& readDistEdges(std::istream& in, ContigGraph<Graph>& g,
		typename graph_traits<Graph>::vertex_descriptor u,
		BetterEP betterEP)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename Graph::edge_property_type EP;
	for (std::string vname; getline(in >> std::ws, vname, ',');) {
		assert(!vname.empty());
		V v = find_vertex(vname, g);
		v = v ^ get(vertex_sense, g, u);
		EP ep;
		in >> ep;
		assert(in);
		if (in.peek() != ' ')
			in >> Ignore(' ');

		E e;
		bool found;
		boost::tie(e, found) = edge(u, v, g);
		if (found) {
			// Parallel edge
			EP& ref = g[e];
			ref = betterEP(ref, ep);
		} else
			g.Graph::add_edge(u, v, ep);
	}
	assert(in.eof());
	return in;
}

/** Read a contig adjacency graph.
 * @param betterEP handle parallel edges
 */
template <typename Graph, typename BetterEP>
std::istream& read_adj(std::istream& in, ContigGraph<Graph>& g,
		BetterEP betterEP)
{
	if (in.eof()) {
		// Allow reading an empty file if the graph is not empty.
		assert(num_vertices(g) > 0);
		return in;
	}
	assert(in);

	typedef typename Graph::vertex_descriptor vertex_descriptor;
	typedef typename Graph::vertex_property_type vertex_property_type;
	typedef typename Graph::edge_property_type edge_property_type;

	// Check for ADJ or DIST format.
	std::string line;
	getline(in, line);
	assert(in);
	unsigned numSemicolons
		= std::count(line.begin(), line.end(), ';');
	if (numSemicolons > 2) {
		std::cerr << "error: expected 0, 1 or 2 semicolons and saw "
			<< numSemicolons << '\n';
		exit(EXIT_FAILURE);
	}
	bool faiFormat = numSemicolons == 0;
	bool adjFormat = numSemicolons == 2;

	// Read the vertex properties.
	if (adjFormat || faiFormat) {
		assert(num_vertices(g) == 0);
		in.clear();
		in.seekg(0, std::ios::beg);
		assert(in);
		for (std::string uname; in >> uname;) {
			vertex_property_type vp;
			if (faiFormat) {
				unsigned length;
				in >> length;
				put(vertex_length, vp, length);
			} else
				in >> vp;
			in >> Ignore('\n');
			assert(in);
			vertex_descriptor u = add_vertex(vp, g);
			put(vertex_name, g, u, uname);
		}
		assert(in.eof());
	}
	assert(num_vertices(g) > 0);
	g_contigNames.lock();

	if (faiFormat)
		return in;

	// Read the edges.
	in.clear();
	in.seekg(0, std::ios::beg);
	assert(in);
	for (std::string name; in >> name;) {
		if (adjFormat)
			in >> Ignore(';');
		vertex_descriptor u = find_vertex(name, false, g);
		for (int sense = false; sense <= true; ++sense) {
			std::string s;
			getline(in, s, !sense ? ';' : '\n');
			assert(in.good());
			std::istringstream ss(s);
			if (!adjFormat) {
				readDistEdges(ss, g, u ^ sense, betterEP);
			} else
			for (std::string vname; ss >> vname;) {
				ss >> std::ws;
				vertex_descriptor v = find_vertex(vname, g);
				assert(!edge(u ^ sense, v ^ sense, g).second);
				if (ss.peek() == '[') {
					ss.get();
					edge_property_type ep;
					ss >> ep >> Ignore(']');
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
