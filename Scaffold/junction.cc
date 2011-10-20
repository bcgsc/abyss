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
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/ref.hpp>
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
using boost::adaptors::filtered;
using boost::cref;
using boost::for_each;
using boost::ref;
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
"  -i, --ignore=FILE     ignore junctions seen in FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by DotIO

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

static inline
void put(vertex_length_t, NoProperty&, unsigned)
{
}

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
	if (opt::verbose > 1)
		printGraphStats(cerr, g);
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
	ContigID::lock();
}

/** Mark contigs that have been seen previously. */
static void markSeen(const string& filePath, vector<bool>& marked)
{
	if (filePath.empty())
		return;
	ifstream in(filePath.c_str());
	assert_good(in, filePath);
	string id;
	ContigPath path;
	while (in >> id >> path) {
		if (path.empty())
			marked[ContigID(id)] = true;
		for (ContigPath::const_iterator it = path.begin();
				it != path.end(); ++it)
			if (!it->ambiguous())
				marked[ContigID(*it)] = true;
	}
	assert(in.eof());
}

/** Return the value of the bit at the specified index. */
struct Marked : unary_function<ContigNode, bool> {
	typedef vector<bool> Data;
	Marked(const Data& data) : m_data(data) { }
	bool operator()(ContigNode u) const
	{
		return m_data[ContigID(u)];
	}
  private:
	const Data& m_data;
};

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
	} else
		readGraph("-", scaffoldG);

	// Add any missing complementary edges.
	addComplementaryEdges(scaffoldG);

	// Read the set of contigs to ignore.
	vector<bool> seen(num_vertices(overlapG) / 2);
	markSeen(opt::ignorePath, seen);
	
	// Extend the junction contigs.
	for_each(vertices(overlapG)
			| filtered(!bind(Marked(seen), _1)),
		bind(extendJunction, cref(overlapG), cref(scaffoldG), _1));

	assert_good(cout, "stdout");

	if (opt::verbose > 0)
		cerr << "Junctions: " << g_count.junctions << '\n'
			<< "Supported: " << g_count.supported << '\n';

	return 0;
}
