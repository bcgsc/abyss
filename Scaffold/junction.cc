#include "config.h"
#include "Common/ContigPath.h"
#include "Common/IOUtil.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include "Uncompress.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/ref.hpp>
#include <algorithm>
#include <cassert>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

using namespace std;
using namespace boost::lambda;
using boost::tie;
#if !__GXX_EXPERIMENTAL_CXX0X__
using boost::cref;
using boost::ref;
#endif

#define PROGRAM "abyss-junction"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2012 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... OVERLAP [SCAFFOLD]...\n"
"Extend junction contigs that are supported by the scaffold graph.\n"
"  OVERLAP   the overlap graph\n"
"  SCAFFOLD  a scaffold graph\n"
"\n"
"  -i, --ignore=FILE     ignore junctions seen in FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by DotIO

	/** Do not check for evidence in the scaffold graph. */
	bool noScaffoldGraph;

	/** Junction contigs to ignore. */
	string ignorePath;

	/** Verbose output. */
	int verbose; // used by PopBubbles

 	/** Output format */
 	int format = DOT; // used by DistanceEst
}

static const char shortopts[] = "i:v";

enum { OPT_HELP = 1, OPT_VERSION, };

static const struct option longopts[] = {
	{ "ignore",      no_argument,       NULL, 'i' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Counts. */
static struct {
	unsigned junctions;
	unsigned supported;
} g_count;

/** An overlap graph. */
typedef DirectedGraph<NoProperty, NoProperty> DG;
typedef ContigGraph<DG> Graph;
typedef Graph OverlapGraph;

/** A scaffold graph. */
typedef Graph ScaffoldGraph;

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
	if (opt::noScaffoldGraph
			|| edge(u, w, scaffoldG).second) {
		// This junction contig is supported by the scaffold graph.
		ContigPath path;
		path.reserve(3);
		extend(overlapG, ~v, back_inserter(path));
		reverseComplement(path.begin(), path.end());
		path.push_back(v);
		extend(overlapG, v, back_inserter(path));
		assert(path.size() >= 3);
		cout << createContigName() << '\t' << path << '\n';
		g_count.supported++;
	}
	g_count.junctions++;
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
	if (opt::verbose > 1)
		printGraphStats(cerr, g);
	g_contigNames.lock();
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'i': arg >> opt::ignorePath; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
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
		for_each(argv + optind, argv + argc,
				bind(readGraph, _1, ref(scaffoldG)));
		// Add any missing complementary edges.
		size_t numAdded = addComplementaryEdges<DG>(scaffoldG);
		if (opt::verbose > 0)
			cerr << "Added " << numAdded << " complementary edges.\n";
		if (opt::verbose > 1)
			printGraphStats(cerr, scaffoldG);
	} else
		opt::noScaffoldGraph = true;

	// Read the set of contigs to ignore.
	vector<bool> seen(num_vertices(overlapG) / 2);
	if (!opt::ignorePath.empty()) {
		ifstream in(opt::ignorePath.c_str());
		assert_good(in, opt::ignorePath);
		markSeenInPath(in, seen);
	}

	// Extend the junction contigs.
	graph_traits<OverlapGraph>::vertex_iterator uit, ulast;
	for (tie(uit, ulast) = vertices(overlapG); uit != ulast; ++uit)
		if (!seen[ContigID(*uit)])
			extendJunction(overlapG, scaffoldG, *uit);

	assert_good(cout, "stdout");

	if (opt::verbose > 0) {
		cerr << "Junctions: " << g_count.junctions << '\n';
		if (!opt::noScaffoldGraph)
			cerr << "Supported: " << g_count.supported << '\n';
	}

	return 0;
}
