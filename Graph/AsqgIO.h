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
	typedef typename vertex_property<Graph>::type VP;
	typedef typename edge_property<Graph>::type EP;

	while (in && in.peek() != EOF) {
		switch (in.peek()) {
		  case 'H':
			in >> expect("HT") >> ignore('\n');
			assert(in);
			break;
		  case 'V': {
			string uname;
			in >> expect("VT") >> uname >> std::ws >> ignore('\t');
			assert(in);
			VP vp;
			put(vertex_length, vp, in.gcount() - 1);
			V u = add_vertex(vp, g);
			put(vertex_name, g, u, uname);
			in >> ignore('\n');
			assert(in);
			break;
		  }
		  case 'E': {
			ContigID u, v;
			unsigned s1, e1, l1, s2, e2, l2, nd;
			bool rc;
			in >> expect("ED") >> u >> v
				>> s1 >> e1 >> l1
				>> s2 >> e2 >> l2
				>> rc >> nd >> ignore('\n');
			assert(in);
			assert(e1 - s1 == e2 - s2);
			assert(((s1 > 0) != (s2 > 0)) ^ rc);
			int d = -(e1 - s1 + 1);
			assert(d < 0);
			add_edge(V(u, s1 == 0), V(v, s2 > 0), EP(d), g);
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
