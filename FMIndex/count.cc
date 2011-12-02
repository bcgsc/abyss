#include "BitUtil.h"
#include "DAWG.h"
#include <algorithm>
#include <boost/graph/depth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <cassert>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

using boost::default_color_type;

#define PROGRAM "abyss-count"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FASTA]\n"
"Count k-mer of the specified file.\n"
"The index file TARGET.fm will be used if present.\n"
"\n"
"  -k, --kmer              the size of a k-mer\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** The size of a k-mer. */
	static unsigned k;

	/** Verbose output. */
	static int verbose;
};

static const char shortopts[] = "k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer", no_argument, NULL, 'k' },
	{ "verbose", no_argument, NULL, 'v' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** A directed acyclic word graph. */
typedef DAWG Graph;

/** A color map that always returns white. */
struct WhiteMap
{
	typedef graph_traits<Graph>::vertex_descriptor key_type;
	typedef default_color_type value_type;
	typedef value_type& reference;
	typedef boost::lvalue_property_map_tag category;
};

default_color_type get(const WhiteMap&,
		graph_traits<Graph>::vertex_descriptor)
{
	return boost::white_color;
}

void put(WhiteMap&,
		graph_traits<Graph>::vertex_descriptor, default_color_type)
{
}

/** Count k-mer. */
class CountKmerVisitor : public boost::default_dfs_visitor
{
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef graph_traits<Graph>::vertex_descriptor V;

  public:
	CountKmerVisitor(vector<char>& s) : m_s(s)
	{
		assert(m_s.empty());
		m_s.reserve(opt::k + 1);
	}

	void discover_vertex(V, const Graph&)
	{
		assert(m_s.size() <= opt::k);
		m_s.push_back('x');
	}

	void finish_vertex(V, const Graph&)
	{
		assert(!m_s.empty());
		m_s.pop_back();
	}

	void examine_edge(E e, const Graph& g)
	{
		m_s.back() = get(boost::edge_name, g, e);
		if (m_s.size() == opt::k) {
			V v = target(e, g);
			unsigned count = v.second - v.first;
			copy(m_s.rbegin(), m_s.rend(),
					ostream_iterator<char>(cout));
			cout << '\t' << count << '\n';
		}
	}

	/** Terminate the search at a depth of k. */
	bool operator()(V, const Graph&) const
	{
		return m_s.size() > opt::k;
	}

  private:
	vector<char>& m_s;
}; // CountKmerVisitor

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
	checkPopcnt();

	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case 'k':
			arg >> opt::k;
			assert(arg.eof());
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
	}

	if (opt::k <= 0) {
		cerr << PROGRAM ": " << "missing -k,--kmer option\n";
		die = true;
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

	// Read the FM index.
	Graph g;
	readFMIndex(g, faPath);

	// Count k-mer.
	vector<char> s;
	boost::depth_first_visit(g, *vertices(g).first,
			CountKmerVisitor(s), WhiteMap(), CountKmerVisitor(s));
	assert_good(cout, "stdout");

	return 0;
}
