#ifndef ASQGIO_H
#define ASQGIO_H 1

#include "Common/ContigID.h"
#include "Common/IOUtil.h"
#include "Graph/Properties.h"
#include <boost/graph/graph_traits.hpp>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

/** Read a graph in ASQG format. */
template <typename Graph>
std::istream& read_asqg(std::istream& in, Graph& g)
{
	assert(in);

	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename vertex_property<Graph>::type VP;
	typedef typename edge_property<Graph>::type EP;

	for (std::string type; in >> type;) {
		if (type == "HT") {
			in >> ignore('\n');
			assert(in);
		} else if (type == "VT") {
			ContigID u;
			std::string seq;
			in >> u >> seq >> ignore('\n');
			assert(in);
			VP vp;
			put(vertex_length, vp, seq.size());
			V x = add_vertex(vp, g);
			assert(u == x);
		} else if (type == "ED") {
			ContigID u, v;
			unsigned s1, e1, l1, s2, e2, l2, nd;
			bool rc;
			in >> u >> v
				>> s1 >> e1 >> l1
				>> s2 >> e2 >> l2
				>> rc >> nd >> ignore('\n');
			assert(in);
			assert(e1 - s1 == e2 - s2);
			assert(((s1 > 0) != (s2 > 0)) ^ rc);
			int d = -(e1 - s1 + 1);
			assert(d < 0);
			add_edge(V(u, s1 == 0), V(v, s2 > 0), EP(d), g);
		} else {
			std::cerr << "error: unknown record type: `"
				<< type << "'\n";
			exit(EXIT_FAILURE);
		}
	}
	assert(in.eof());
	return in;
}

#endif
