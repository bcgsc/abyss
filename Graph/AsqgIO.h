#ifndef ASQGIO_H
#define ASQGIO_H 1

#include "Common/IOUtil.h"
#include "Graph/Properties.h"
#include <boost/graph/graph_traits.hpp>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

using boost::graph_traits;

/** Write a graph in ASQG format. */
template <typename Graph>
std::ostream& write_asqg(std::ostream& out, Graph& g)
{
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename graph_traits<Graph>::edge_iterator Eit;
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::vertex_iterator Vit;
	typedef typename vertex_bundle_type<Graph>::type VP;

	out << "HT\tVN:i:1\n";
	assert(out);

	std::pair<Vit, Vit> vrange = vertices(g);
	for (Vit uit = vrange.first; uit != vrange.second; ++uit, ++uit) {
		V u = *uit;
		if (get(vertex_removed, g, u))
			continue;
		const VP& vp = g[u];
		out << "VT\t" << get(vertex_contig_name, g, u)
			<< "\t*\tLN:i:" << vp.length;
		if (vp.coverage > 0)
			out << "\tKC:i:" << vp.coverage;
		out << '\n';
	}

	std::pair<Eit, Eit> erange = edges(g);
	for (Eit eit = erange.first; eit != erange.second; ++eit) {
		E e = *eit;
		V u = source(e, g);
		V v = target(e, g);
		if (get(vertex_removed, g, u))
			continue;

		// Output only the canonical edge.
		if (u > get(vertex_complement, g, v))
			continue;

		assert(!get(vertex_removed, g, v));
		int distance = g[e].distance;
		assert(distance <= 0);
		unsigned overlap = -distance;
		unsigned ulen = g[u].length;
		unsigned vlen = g[v].length;
		bool usense = get(vertex_sense, g, u);
		bool vsense = get(vertex_sense, g, v);
		out << "ED\t" << get(vertex_contig_name, g, u)
			<< ' ' << get(vertex_contig_name, g, v)
			<< ' ' << (usense ? 0 : ulen - overlap)
			<< ' ' << int((usense ? overlap : ulen) - 1)
			<< ' ' << ulen
			<< ' ' << (!vsense ? 0 : vlen - overlap)
			<< ' ' << int((!vsense ? overlap : vlen) - 1)
			<< ' ' << vlen
			<< ' ' << (usense != vsense)
			<< " -1\n"; // number of mismatches
	}
	return out;
}

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
			in >> expect("VT") >> uname >> seq;
			assert(in);
			assert(!seq.empty());

			unsigned length = 0;
			if (seq == "*") {
				in >> expect(" LN:i:") >> length;
				assert(in);
			} else
				length = seq.size();

			unsigned coverage = 0;
			if (in.peek() == '\t' && in.get() == '\t' && in.peek() == 'K') {
				in >> expect("KC:i:") >> coverage;
				assert(in);
			}

			in >> Ignore('\n');
			assert(in);

			if (addVertices) {
				VP vp;
				put(vertex_length, vp, length);
				put(vertex_coverage, vp, coverage);
				V u = add_vertex(vp, g);
				put(vertex_name, g, u, uname);
			} else {
				V u = find_vertex(uname, false, g);
				assert(get(vertex_index, g, u) < num_vertices(g));
				(void)u;
			}
			break;
		  }
		  case 'E': {
			std::string uname, vname;
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
			V u = find_vertex(uname, s1 == 0, g);
			V v = find_vertex(vname, s2 > 0, g);
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
