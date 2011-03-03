#include "config.h"
#include "ContigGraph.h"
#include "ContigGraphAlgorithms.h"
#include "ContigNode.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "DirectedGraph.h"
#include "DotIO.h"
#include "GraphUtil.h"
#include "IOUtil.h"
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

using namespace std;

#define PROGRAM "abyss-scaffold"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [DIST]\n"
"Scaffold contigs using the distance estimate graph.\n"
"  DIST  estimates of the distance between contigs\n"
"\n"
"  -n, --npairs=N        minimum number of pairs [0]\n"
"  -s, --seed-length=N   minimum contig length [0]\n"
"  -o, --out=FILE        write the paths to FILE\n"
"  -g, --graph=FILE      write the graph to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties
	int dot; // used by Estimate

	/** Minimum number of pairs. */
	static unsigned minNumPairs;

	/** Minimum contig length. */
	static unsigned minContigLength;

	/** Write the paths to this file. */
	static string out;

	/** Write the graph to this file. */
	static string graphPath;

	/** Verbose output. */
	static int verbose;
}

static const char shortopts[] = "g:n:o:s:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MAX_COST };

static const struct option longopts[] = {
	{ "graph",       no_argument,       NULL, 'g' },
	{ "npairs",      required_argument, NULL, 'n' },
	{ "out",         required_argument, NULL, 'o' },
	{ "seed-length", required_argument, NULL, 's' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Contig length property. */
struct Length {
	unsigned length;

	Length& operator+=(const Length& o)
	{
		length += o.length;
		return *this;
	}

	Length& operator+=(const Distance& o)
	{
		assert((int)length + (int)o.distance > 0);
		length += o.distance;
		return *this;
	}

	friend ostream& operator<<(ostream& out, const Length& o)
	{
		return out << "l=" << o.length;
	}

	friend istream& operator>>(istream& in, Length& o)
	{
		return in >> expect("l =") >> o.length;
	}
};

/** A distance estimate graph. */
typedef DirectedGraph<Length, Distance> DG;
typedef ContigGraph<DG> Graph;

/** Add missing complementary edges. */
static void addComplementaryEdges(DG& g)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::edge_descriptor E;
	typedef GTraits::edge_iterator Eit;
	typedef GTraits::vertex_descriptor V;

	std::pair<Eit, Eit> erange = edges(g);
	unsigned numAdded;
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

/** Remove short vertices and unsupported edges from the graph. */
static void filterGraph(Graph& g)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_descriptor V;
	typedef GTraits::vertex_iterator Vit;

	unsigned numRemovedV = 0;
	std::pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		if (g[u].length < opt::minContigLength)
			clear_vertex(u, g);
		if (out_degree(u, g) == 0 && in_degree(u, g) == 0) {
			remove_vertex(u, g);
			numRemovedV++;
		}
	}
	if (opt::verbose > 0)
		cerr << "Removed " << numRemovedV << " vertices.\n";
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'g': arg >> opt::graphPath; break;
			case 'n': arg >> opt::minNumPairs; break;
			case 'o': arg >> opt::out; break;
			case 's': arg >> opt::minContigLength; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (argc - optind < 0) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	} else if (argc - optind > 1) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	string distFilePath(optind < argc ? argv[optind++] : "-");

	// Read the distance estimate graph.
	ifstream fin(distFilePath.c_str());
	istream& in = distFilePath == "-" ? cin : fin;
	assert_good(in, distFilePath);
	Graph g;
	read_dot<DG>(in, g);
	assert(in.eof());
	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	// Filter the graph.
	addComplementaryEdges(g);
	filterGraph(g);

	// Output the graph.
	if (!opt::graphPath.empty()) {
		ofstream out(opt::graphPath.c_str());
		assert_good(out, opt::graphPath);
		write_dot(out, g);
		assert_good(out, opt::graphPath);
	}

	// Assemble the paths.
	typedef vector<ContigPath> ContigPaths;
	ContigPaths paths;
	assemble(g, back_inserter(paths));
	sort(paths.begin(), paths.end());

	// Output the paths.
	ofstream fout(opt::out.c_str());
	ostream& out = opt::out.empty() || opt::out == "-" ? cout : fout;
	assert_good(out, opt::out);
	for (vector<ContigPath>::const_iterator it = paths.begin();
			it != paths.end(); ++it)
		out << ContigID::create() << '\t' << *it << '\n';
	assert_good(out, opt::out);

	return 0;
}
