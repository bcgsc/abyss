/** Convert a graph from adj format to dot format.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 * Copyright 2010 Genome Sciences Centre
 */
#include "ContigGraph.h"
#include "DirectedGraph.h"
#include <fstream>
#include <getopt.h>
#include <iostream>

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
"      --help     display this help and exit\n"
"      --version  output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static int k;
}

static const char shortopts[] = "k:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",    required_argument, NULL, 'k' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Contig properties. */
struct Contig {
	unsigned length;
	unsigned coverage;

	Contig() { }
	Contig(unsigned length, unsigned coverage)
		: length(length), coverage(coverage) { }

	friend ostream& operator <<(ostream& out, const Contig& o)
	{
		out << "l=" << o.length;
		if (opt::k > 0)
			return out << " c="
				<< (float)o.coverage / (o.length - opt::k + 1);
		else
			return out << " C=" << o.coverage;
	}

	friend istream& operator >>(istream& in, Contig& o)
	{
		return in >> o.length >> o.coverage;
	}
};

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?': die = true; break;
		  case 'k': arg >> opt::k; break;
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

	typedef ContigGraph<DirectedGraph<Contig> > Graph;
	Graph g;
	in >> g;
	assert(in.eof());

	cout << "digraph \"" << path << "\" {\n";
	if (opt::k > 0)
		cout << "k=" << opt::k << "\n"
			"edge[d=" << -(opt::k-1) << "]\n";
	cout << dot_writer<Graph>(g) << "}\n";
	assert(cout.good());
	return 0;
}
