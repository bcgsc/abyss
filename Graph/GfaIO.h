#ifndef GFAIO_H
#define GFAIO_H 1

#include "Common/IOUtil.h"
#include "Graph/Properties.h"
#include <boost/graph/graph_traits.hpp>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

using boost::graph_traits;

/** Write a graph in GFA format. */
template <typename Graph>
std::ostream& write_gfa(std::ostream& out, Graph& g)
{
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename graph_traits<Graph>::edge_iterator Eit;
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::vertex_iterator Vit;
	typedef typename vertex_bundle_type<Graph>::type VP;

	out << "H\tVN:Z:1.0\n";
	assert(out);

	std::pair<Vit, Vit> vrange = vertices(g);
	for (Vit uit = vrange.first; uit != vrange.second; ++uit, ++uit) {
		V u = *uit;
		if (get(vertex_removed, g, u))
			continue;
		const VP& vp = g[u];
		out << "S\t" << get(vertex_contig_name, g, u)
			<< "\t*\t*\tLN:i:" << vp.length;
		if (vp.coverage > 0)
			out << "\tXC:i:" << vp.coverage;
		out << '\n';
	}

	std::pair<Eit, Eit> erange = edges(g);
	for (Eit eit = erange.first; eit != erange.second; ++eit) {
		E e = *eit;
		V u = source(e, g);
		V v = target(e, g);
		if (v < u || get(vertex_removed, g, u))
			continue;
		assert(!get(vertex_removed, g, v));
		int distance = g[e].distance;
		out << 'L'
			<< '\t' << get(vertex_name, g, u)
			<< '\t' << get(vertex_name, g, v);
		if (distance < 0)
			out << '\t' << -distance << "M\n";
		else
			out << "\t*\n";
		assert(out);
	}
	return out;
}

/** Read a graph in GFA format. */
template <typename Graph>
std::istream& read_gfa(std::istream& in, Graph& g)
{
	assert(in);

	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename vertex_property<Graph>::type VP;
	typedef typename edge_property<Graph>::type EP;

	// Add vertices if this graph is empty.
	bool addVertices = num_vertices(g) == 0;

	while (in && in.peek() != EOF) {
		switch (in.peek()) {
		  case 'H':
			in >> expect("H\tVN:Z:1.0\n");
			assert(in);
			break;

		  case 'S': {
			std::string uname, seq;
			in >> expect("S\t") >> uname >> seq;
			assert(in);
			assert(!seq.empty());

			unsigned length = 0;
			if (seq == "*") {
				in >> expect(" LN:i:") >> length;
				assert(in);
			} else
				length = seq.size();
			in >> Ignore('\n');
			assert(in);

			if (addVertices) {
				VP vp;
				put(vertex_length, vp, length);
				V u = add_vertex(vp, g);
				put(vertex_name, g, u, uname);
			} else {
				V u = find_vertex(uname, false, g);
				assert(get(vertex_index, g, u) < num_vertices(g));
				(void)u;
			}
			break;
		  }

		  case 'L': {
			std::string uname, vname;
			unsigned overlap;
			in >> expect("L\t") >> uname >> vname >> std::ws;
			if (in.peek() == '*') {
				in.get();
				overlap = 0;
			} else {
				in >> overlap >> expect("M");
			}
			in >> Ignore('\n');
			assert(in);
			V u = find_vertex(uname, g);
			V v = find_vertex(vname, g);
			if (overlap > 0) {
				int d = -overlap;
				EP ep(d);
				add_edge(u, v, ep, g);
			} else
				add_edge(u, v, g);
			break;
		  }

		  default: {
			std::string s;
			in >> s;
			std::cerr << "error: unknown record type: `"
				<< s << "'\n";
			exit(EXIT_FAILURE);
		  }
		}
	}
	assert(in.eof());
	return in;
}

#endif
