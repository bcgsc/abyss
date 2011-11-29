#ifndef GRAPHIO_H
#define GRAPHIO_H 1

#include "Graph/Options.h"
#include "AdjIO.h"
#include "AsqgIO.h"
#include "DistIO.h"
#include "DotIO.h"
#include "FastaIO.h"
#include "SAMIO.h"
#include <cassert>
#include <cstdlib> // for abort
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
	  case ASQG:
		return write_asqg(out, g);
	  case DIST:
		return write_dist(out, g);
	  case DOT: case DOT_MEANCOV:
		return out << dot_writer(g);
	  case SAM:
		return write_sam(out, g, program, commandLine);
	  default:
		assert(false);
		abort();
	}
}

#include "ContigGraph.h"

/** Read a graph. */
template <typename Graph, typename BetterEP>
std::istream& read_graph(std::istream& in, ContigGraph<Graph>& g,
		BetterEP betterEP)
{
	in >> std::ws;
	assert(in);
	switch (in.peek()) {
	  case 'd': // digraph: GraphViz dot format
		return read_dot<Graph>(in, g, betterEP);
	  case 'H': // HT: ASQG format
		return read_asqg(in, g);
	  case '>': // FASTA format for vertices
		return read_fasta(in, g);
	  default: // adj format
		return read_adj(in, g, betterEP);
	}
}

/** Disallow parallel edges. */
struct DisallowParallelEdges {
	template <typename EP>
	EP operator()(const EP& a, const EP& b) const
	{
		std::cerr << "error: parallel edges:"
			" [" << a << "], [" << b << "]\n";
		exit(EXIT_FAILURE);
	}
};

/** Read a graph. */
template <typename Graph>
std::istream& operator>>(std::istream& in, ContigGraph<Graph>& g)
{
	return read_graph(in, g, DisallowParallelEdges());
}

#endif
