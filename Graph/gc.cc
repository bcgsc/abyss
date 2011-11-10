/** Convert a graph to dot format.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */
#include "ContigGraph.h"
#include "DirectedGraph.h"
#include "GraphIO.h"
#include "GraphUtil.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include <fstream>
#include <getopt.h>
#include <iostream>

using namespace std;

#define PROGRAM "abyss-gc"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [FILE]...\n"
"Count the number of vertices and edges in a graph.\n"
"\n"
"  -v, --verbose  display verbose output\n"
"      --help     display this help and exit\n"
"      --version  output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by DotIO
	static int verbose;
}

static const char shortopts[] = "v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose", no_argument, NULL, 'v' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** No property. */
struct NoProperty
{
	NoProperty(...) { }
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

template <typename Tag>
void put(Tag, NoProperty&, unsigned)
{
}

/** Read a graph from the specified file. */
template <typename Graph, typename BetterEP>
static void readGraph(const string& path, Graph& g, BetterEP betterEP)
{
	if (opt::verbose > 0)
		cout << "Reading `" << path << "'...\n";
	ifstream fin(path.c_str());
	istream& in = path == "-" ? cin : fin;
	assert_good(in, path);
	read_graph(in, g, betterEP);
	assert(in.eof());
	printGraphStats(cout, g);
	ContigID::lock();
}

/** Read a graph from the specified files. */
template <typename Graph, typename It, typename BetterEP>
void readGraphs(Graph& g, It first, It last, BetterEP betterEP)
{
	if (first != last) {
		for (It it = first; it < last; ++it)
			readGraph(*it, g, betterEP);
	} else
		readGraph("-", g, betterEP);
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?': die = true; break;
		  case 'v': opt::verbose++; break;
		  case OPT_HELP:
			cout << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case OPT_VERSION:
			cout << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	ContigGraph<DirectedGraph<NoProperty, NoProperty> > g;
	readGraphs(g, argv + optind, argv + argc,
			DisallowParallelEdges());
	assert_good(cout, "-");

	return 0;
}
