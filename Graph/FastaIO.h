#ifndef FASTAIO_H
#define FASTAIO_H 1

#include "Common/ContigID.h"
#include "Graph/Properties.h"
#include <boost/graph/graph_traits.hpp>
#include <cassert>
#include <istream>
#include <sstream>
#include <string>

/** Read the vertices from a FASTA file. */
template <typename Graph>
std::istream& read_fasta(std::istream& in, Graph& g)
{
	assert(in);

	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename vertex_property<Graph>::type VP;

	for (ContigID u; in.peek() != EOF && in >> expect(">") >> u;) {
		std::string comment;
		getline(in, comment);
		assert(in);
		std::istringstream ss(comment);
		VP vp;
		ss >> vp;
		size_t n = 0;
		while (in >> std::ws && in.peek() != '>' && in) {
			in >> ignore('\n');
			assert(in);
			assert(in.gcount() > 1);
			n += in.gcount() - 1;
		}
		put(vertex_length, vp, n);
		V x = add_vertex(vp, g);
		assert(u == x);
	}
	assert(in.eof());
	return in;
}

#endif
