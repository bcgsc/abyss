/** Convert a graph from adj format to dot format.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 * Copyright 2010 Genome Sciences Centre
 */
#include "ContigGraph.h"
#include "ContigProperties.h"
#include "DirectedGraph.h"
#include "GraphIO.h"
#include "Histogram.h"
#include "Uncompress.h"
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator> // for ostream_iterator

using namespace std;

#define PROGRAM "abyss-adjtodot"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " FILE\n"
"Convert the specified graph from adj format to dot format.\n"
"\n"
"  -k, --kmer=N   report the mean k-mer coverage, otherwise\n"
"                 the sum k-mer coverage is reported\n"
"      --adj             output the results in adj format\n"
"      --dot             output the results in dot format [default]\n"
"      --sam             output the results in SAM format\n"
"  -v, --verbose  display verbose output\n"
"      --help     display this help and exit\n"
"      --version  output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
 	unsigned k; // used by Distance
	static int verbose;

	/** Output format */
	int format = DOT; // used by ContigProperties
}

static const char shortopts[] = "k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "adj",     no_argument,       &opt::format, ADJ },
	{ "dot",     no_argument,       &opt::format, DOT },
	{ "sam",     no_argument,       &opt::format, SAM },
	{ "kmer",    required_argument, NULL, 'k' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

int main(int argc, char** argv)
{
	string commandLine;
	{
		ostringstream ss;
		char** last = argv + argc - 1;
		copy(argv, last, ostream_iterator<const char *>(ss, " "));
		ss << *last;
		commandLine = ss.str();
	}

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

	string path = optind < argc ? argv[optind] : "-";

	ifstream fin(path.c_str());
	if (path != "-")
		assert(fin.is_open());
	istream& in = path == "-" ? cin : fin;

	typedef ContigGraph<DirectedGraph<
		ContigProperties, Distance> > Graph;
	Graph g;
	in >> g;
	assert(in.eof());

	if (opt::verbose > 0) {
		unsigned v = num_vertices(g);
		unsigned e = num_edges(g);
		cerr << "V=" << v
			<< " E=" << e
			<< " E/V=" << (float)e / v
			<< endl;

		// Print a histogram of the out degree of vertices.
		Histogram h;
		typedef Graph::vertex_iterator vertex_iterator;
		std::pair<vertex_iterator, vertex_iterator>
			vit = vertices(g);
		for (vertex_iterator u = vit.first; u != vit.second; ++u)
			h.insert(out_degree(*u, g));
		cerr <<
			"Degree: " << h.barplot() << "\n"
			"        01234" << endl;
	}

	write_graph(cout, g, PROGRAM, commandLine);
	assert(cout.good());
	return 0;
}
