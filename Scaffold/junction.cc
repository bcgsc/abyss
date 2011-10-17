#include "config.h"
#include "Common/ContigPath.h"
#include "Common/IOUtil.h"
#include "Graph/ContigGraph.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <algorithm>
#include <cassert>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

using namespace std;
using namespace boost::lambda;
using boost::range::for_each;
using boost::tie;

#define PROGRAM "abyss-junction"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... OVERLAP [SCAFFOLD]...\n"
"Extend junction contigs that are supported by the scaffold graph.\n"
"  OVERLAP   the overlap graph\n"
"  SCAFFOLD  a scaffold graph\n"
"\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by DotIO

	/** Verbose output. */
	int verbose; // used by PopBubbles

 	/** Output format */
 	int format = DOT; // used by DistanceEst
}

static const char shortopts[] = "v";

enum { OPT_HELP = 1, OPT_VERSION, };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** No property. */
struct NoProperty {
	bool operator==(const NoProperty&) const { return true; }
	friend ostream& operator<<(ostream& out, const NoProperty&)
	{
		return out;
	}
	friend istream& operator>>(istream& in, NoProperty&)
	{
		return in;
	}
};

/** An overlap graph. */
typedef DirectedGraph<NoProperty, NoProperty> DG;
typedef ContigGraph<DG> Graph;
typedef Graph OverlapGraph;

/** A scaffold graph. */
typedef Graph ScaffoldGraph;

/** Add missing complementary edges. */
static void addComplementaryEdges(DG& g)
{
	typedef graph_traits<DG> GTraits;
	typedef GTraits::edge_descriptor E;
	typedef GTraits::edge_iterator Eit;
	typedef GTraits::vertex_descriptor V;

	pair<Eit, Eit> erange = edges(g);
	unsigned numAdded = 0;
	for (Eit eit = erange.first; eit != erange.second; ++eit) {
		E e = *eit;
		V u = source(e, g), v = target(e, g);
		if (!edge(~v, ~u, g).second) {
			add_edge(~v, ~u, g[e], g);
			numAdded++;
		}
	}
	if (opt::verbose > 0)
		cerr << "Added " << numAdded << " complementary edges.\n";
}

/** Extend junction contigs. */
void extendJunction(
		const OverlapGraph& overlapG,
		const ScaffoldGraph& scaffoldG,
		graph_traits<OverlapGraph>::vertex_descriptor v)
{
	if (get(vertex_sense, overlapG, v)
			|| in_degree(v, overlapG) != 1
			|| out_degree(v, overlapG) != 1)
		return;
	typedef graph_traits<OverlapGraph>::vertex_descriptor V;
	V u = source(*in_edges(v, overlapG).first, overlapG);
	V w = *adjacent_vertices(v, overlapG).first;
	if (edge(u, v, scaffoldG).second
			&& edge(v, w, scaffoldG).second
			&& edge(u, w, scaffoldG).second) {
		// This junction contig is supported by the scaffold graph.
		cout << ContigID::create() << '\t'
			<< u << ' ' << v << ' ' << w << '\n';
	}
}

/** Allow parallel edges. */
struct AllowParallelEdges {
	template <typename EP>
	EP operator()(const EP& a, const EP&) const
	{
		return a;
	}
};

/** Read a graph from the specified file. */
static void readGraph(const string& path, Graph& g)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << path << "'...\n";
	ifstream fin(path.c_str());
	istream& in = path == "-" ? cin : fin;
	assert_good(in, path);
	read_graph(in, g, AllowParallelEdges());
	assert(in.eof());
	if (opt::verbose > 0)
		printGraphStats(cerr, g);
	ContigID::lock();
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	OverlapGraph overlapG;
	readGraph(argv[optind++], overlapG);

	ScaffoldGraph scaffoldG(overlapG.num_vertices() / 2);
	if (optind < argc) {
		for (; optind < argc; optind++)
			readGraph(argv[optind], scaffoldG);
	} else
		readGraph("-", scaffoldG);

	// Add any missing complementary edges.
	addComplementaryEdges(scaffoldG);

	for_each(vertices(overlapG),
		bind(extendJunction, overlapG, scaffoldG, _1));

	assert_good(cout, "stdout");

	return 0;
}
