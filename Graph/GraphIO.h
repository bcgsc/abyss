#ifndef GRAPHIO_H
#define GRAPHIO_H 1

#include "Graph/Options.h"
#include "AdjIO.h"
#include "DotIO.h"
#include "SAMIO.h"
#include <cassert>
#include <istream>
#include <ostream>
#include <string>

/** Write the graph g to the specified output stream in the format
 * specified by opt::format.
 */
template <typename Graph>
std::ostream& write_graph(std::ostream& out, const Graph& g,
		const std::string& program, const std::string& commandLine)
{
	switch (opt::format) {
	  case ADJ:
		return out << adj_writer<Graph>(g);
	  case DOT:
		return out << dot_writer(g);
	  case SAM:
		out << "@HD\tVN:1.0\n"
			"@PG\tID:" << program << "\tVN:" VERSION "\t"
			"CL:" << commandLine << '\n';
		return out << sam_writer<Graph>(g);
	  default:
		assert(false);
	}
}

#include "ContigGraph.h"

template <typename Graph>
std::istream& operator>>(std::istream& in, ContigGraph<Graph>& g)
{
	in >> std::ws;
	assert(in);
	switch (in.peek()) {
	  case 'd': // digraph: GraphViz dot format
		return read_dot<Graph>(in, g);
	  default: // adj format
		return read_adj(in, g);
	}
}

#endif
