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

using boost::graph_traits;

/** Read a graph in ASQG format. */
template <typename Graph>
std::istream& read_asqg(std::istream& in, Graph& g)
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
			in >> expect("HT") >> Ignore('\n');
			assert(in);
			break;
		  case 'V': {
			std::string uname, seq;
			in >> expect("VT") >> uname >> seq >> Ignore('\n');
			assert(in);
			assert(!seq.empty());

			if (addVertices) {
				VP vp;
				put(vertex_length, vp, seq.size());
				V u = add_vertex(vp, g);
				put(vertex_name, g, u, uname);
			} else {
				V u(uname, false);
				assert(get(vertex_index, g, u) < num_vertices(g));
			}
			break;
		  }
		  case 'E': {
			ContigID uname, vname;
			unsigned s1, e1, l1, s2, e2, l2;
			bool rc;
			int nd;
			in >> expect("ED") >> uname >> vname
				>> s1 >> e1 >> l1
				>> s2 >> e2 >> l2
				>> rc >> nd >> Ignore('\n');
			assert(in);
			assert(s1 < e1 && e1 < l1 && s2 < e2 && e2 < l2);
			assert(e1 - s1 == e2 - s2);
			assert(e1 - s1 + 1 < l1 && e2 - s2 + 1 < l2);
			assert(((s1 > 0) == (s2 > 0)) == rc);
			int d = -(e1 - s1 + 1);
			assert(d < 0);
			EP ep(d);
			V u(uname, s1 == 0), v(vname, s2 > 0);
			std::pair<E, bool> e = edge(u, v, g);
			if (e.second) {
				// Ignore duplicate edges that are self loops.
				assert(g[e.first] == ep);
				assert(u == v);
			} else
				add_edge(u, v, ep, g);
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
