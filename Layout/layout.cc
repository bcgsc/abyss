#include "config.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "Uncompress.h"
#include "Graph/Assemble.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphAlgorithms.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

using namespace std;
using boost::tie;

#define PROGRAM "abyss-layout"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2012 Shaun Jackman\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... OVERLAP\n"
"Layout contigs using the sequence overlap graph.\n"
"Output sequence paths.\n"
"  OVERLAP  the sequence overlap graph\n"
"\n"
"  -s, --min-length=N    minimum sequence length [0]\n"
"  -m, --min-overlap=N   minimum overlap [0]\n"
"  -k, --kmer=N          length of a k-mer\n"
"  -o, --out=FILE        write the paths to FILE\n"
"  -g, --graph=FILE      write the graph to FILE\n"
"      --tred            remove transitive edges\n"
"      --no-tred         do not remove transitive edges [default]\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties

	/** Minimum sequence length. */
	static unsigned minLength;

	/** Minimum overlap. */
	static unsigned minOverlap;

	/** Write the paths to this file. */
	static string out;

	/** Write the graph to this file. */
	static string graphPath;

	/** Remove transitive edges. */
	static int tred;

	/** Verbose output. */
	int verbose; // used by PopBubbles

	/** Output format */
	int format = DOT;
}

static const char shortopts[] = "g:k:m:o:s:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "graph",       required_argument, NULL, 'g' },
	{ "kmer",        required_argument, NULL, 'k' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "out",         required_argument, NULL, 'o' },
	{ "min-length",  required_argument, NULL, 's' },
	{ "tred",        no_argument, &opt::tred, true },
	{ "no-tred",     no_argument, &opt::tred, false },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** An overlap graph. */
typedef DirectedGraph<Length, Distance> DG;
typedef ContigGraph<DG> Graph;

/** Remove short vertices. */
static void filterVertices(Graph& g, unsigned minLength)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_descriptor V;
	typedef GTraits::vertex_iterator Vit;

	if (minLength == 0)
		return;

	// Remove short sequences.
	unsigned numRemoved = 0;
	std::pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		if (g[u].length < minLength)
			clear_vertex(u, g);
		if (out_degree(u, g) == 0 && in_degree(u, g) == 0) {
			remove_vertex(u, g);
			numRemoved++;
		}
	}

	if (opt::verbose > 0) {
		cerr << "Ignored " << numRemoved << " sequences shorter than "
			<< minLength << " bp.\n";
		printGraphStats(cerr, g);
	}
}

/** Return true if the edge is a small overlap. */
struct IsSmallOverlap {
	IsSmallOverlap(Graph& g) : m_g(g) { }
	bool operator()(graph_traits<Graph>::edge_descriptor e) const
	{
		int maxDistance = -opt::minOverlap;
		return m_g[e].distance > maxDistance;
	}
	const Graph& m_g;
};

/** Remove small overlaps. */
static void filterEdges(Graph& g, unsigned minOverlap)
{
	if (minOverlap == 0)
		return;
	unsigned numBefore = num_edges(g);
	remove_edge_if(IsSmallOverlap(g), static_cast<DG&>(g));
	unsigned numRemoved = numBefore - num_edges(g);
	if (opt::verbose > 0) {
		cerr << "Removed " << numRemoved << " small overlaps.\n";
		printGraphStats(cerr, g);
	}
}

/** Read a graph from the specified file. */
static void readGraph(const string& path, Graph& g)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << path << "'...\n";
	ifstream fin(path.c_str());
	istream& in = path == "-" ? cin : fin;
	assert_good(in, path);
	in >> g;
	assert(in.eof());
	if (opt::verbose > 0)
		printGraphStats(cerr, g);
	g_contigNames.lock();
}

/** Return the length histogram. */
static Histogram buildLengthHistogram(const Graph& g)
{
	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef graph_traits<Graph>::vertex_iterator Vit;
	Histogram h;
	Vit uit, ulast;
	for (tie(uit, ulast) = vertices(g); uit != ulast; ++++uit) {
		V u = *uit;
		if (!get(vertex_removed, g, u))
			h.insert(g[u].length);
	}
	return h;
}

/** Run abyss-layout. */
int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			die = true;
			break;
		  case 'k':
			arg >> opt::k;
			break;
		  case 'g':
			arg >> opt::graphPath;
			break;
		  case 'm':
			arg >> opt::minOverlap;
			break;
		  case 'o':
			arg >> opt::out;
			break;
		  case 's':
			arg >> opt::minLength;
			break;
		  case 'v':
			opt::verbose++;
			break;
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

	Graph g;
	if (optind < argc) {
		for (; optind < argc; optind++)
			readGraph(argv[optind], g);
	} else
		readGraph("-", g);

	// Remove short sequences.
	filterVertices(g, opt::minLength);

	// Remove small overlaps.
	filterEdges(g, opt::minOverlap);

	// Remove transitive edges.
	if (opt::tred) {
		unsigned numTransitive = remove_transitive_edges(g);
		if (opt::verbose > 0) {
			cerr << "Removed " << numTransitive
				<< " transitive edges.\n";
			printGraphStats(cerr, g);
		}
	}

	/** A container of contig paths. */
	typedef vector<ContigPath> ContigPaths;

	// Assemble the paths.
	ContigPaths paths;
	assembleDFS(g, back_inserter(paths));
	sort(paths.begin(), paths.end());
	if (opt::verbose > 0) {
		unsigned n = 0;
		for (ContigPaths::const_iterator it = paths.begin();
				it != paths.end(); ++it)
			n += it->size();
		cerr << "Assembled " << n << " sequences in "
			<< paths.size() << " contigs.\n";
		printGraphStats(cerr, g);
	}

	// Output the paths.
	ofstream fout(opt::out.c_str());
	ostream& out = opt::out.empty() || opt::out == "-" ? cout : fout;
	assert_good(out, opt::out);
	g_contigNames.unlock();
	for (vector<ContigPath>::const_iterator it = paths.begin();
			it != paths.end(); ++it)
		out << createContigName() << '\t' << *it << '\n';
	assert_good(out, opt::out);

	// Create the new vertices.
	for (vector<ContigPath>::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		const ContigPath& path = *it;
		merge(g, path.begin(), path.end());
		remove_vertex_if(g, path.begin(), path.end(),
				not1(std::mem_fun_ref(&ContigNode::ambiguous)));
	}
	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	// Output the graph.
	if (!opt::graphPath.empty()) {
		ofstream out(opt::graphPath.c_str());
		assert_good(out, opt::graphPath);
		write_dot(out, g);
		assert_good(out, opt::graphPath);
	}

	// Print assembly contiguity statistics.
	if (opt::verbose > 0) {
		Histogram h = buildLengthHistogram(g);
		const unsigned STATS_MIN_LENGTH = 200; // bp
		printContiguityStats(cerr, h, STATS_MIN_LENGTH) << '\n';
	}

	return 0;
}
