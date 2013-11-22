#include "BitUtil.h"
#include "DAWG.h"
#include "Uncompress.h"
#include <algorithm>
#include <boost/graph/depth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <cassert>
#include <cctype>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <map>
#include <vector>

using namespace std;

using boost::default_dfs_visitor;

#define PROGRAM "abyss-dawg"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FASTA]\n"
"Output a directed acyclic word graph (DAWG) of the specified file.\n"
"The index file TARGET.fm will be used if present.\n"
"\n"
" Options:\n"
"\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** Verbose output. */
	static int verbose;
};

static const char shortopts[] = "v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose", no_argument, NULL, 'v' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** A directed acyclic word graph. */
typedef DAWG Graph;

/** DAWG visitor. */
struct DAWGVisitor : public default_dfs_visitor
{
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef graph_traits<Graph>::vertex_descriptor V;

	void examine_edge(E e, const Graph& g)
	{
		using boost::edge_name;
		V u = source(e, g);
		V v = target(e, g);
		char c = get(edge_name, g, e);
		cout << '"'
			<< u.first << ',' << u.second << "\" -> \""
			<< v.first << ',' << v.second << "\""
			" [label=\"";
		if (isprint(c))
			cout << c;
		else
			cout << "0x" << hex << (unsigned)c << dec;
		cout << "\"]\n";
	}
};

/** Read an FM index. */
static void readFMIndex(FMIndex& g, const string& faPath)
{
	string fmPath = faPath + ".fm";
	ifstream in(fmPath.c_str());
	if (in) {
		if (opt::verbose > 0)
			cerr << "Reading `" << fmPath << "'...\n";
		assert_good(in, fmPath);
		in >> g;
		assert_good(in, fmPath);
		in.close();
		return;
	}

	// Build the FM index.
	std::vector<FMIndex::value_type> s;
	if (string(faPath) == "-") {
		if (opt::verbose > 0)
			std::cerr << "Reading stdin...\n";
		copy(istreambuf_iterator<char>(cin),
				istreambuf_iterator<char>(),
				back_inserter(s));
		assert_good(cin, "stdin");
	} else {
		if (opt::verbose > 0)
			std::cerr << "Reading `" << faPath << "'...\n";
		readFile(faPath.c_str(), s);
	}
	g.assign(s.begin(), s.end());
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
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

	if (argc - optind > 1) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}
	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	string faPath(optind < argc ? argv[optind] : "-");

	using boost::default_color_type;
	using boost::depth_first_visit;
	using boost::make_assoc_property_map;

	typedef graph_traits<Graph>::vertex_descriptor V;

	// Read the FM index.
	Graph g;
	readFMIndex(g, faPath);

	cout << "digraph dawg {\n";

	map<V, default_color_type> colorMap;
	depth_first_visit(g, *vertices(g).first,
			DAWGVisitor(), make_assoc_property_map(colorMap));

	cout << "}\n" << flush;
	assert_good(cout, "stdout");

	return 0;
}
