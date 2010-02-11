/**
 * Identify and pop simple bubbles.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "config.h"
#include "Common/Options.h"
#include "AffixIterator.h"
#include "ContigNode.h"
#include "Sequence.h"
#include "Sense.h"
#include <algorithm>
#include <cerrno>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits> // for numeric_limits
#include <map>
#include <set>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

#define PROGRAM "PopBubbles"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... ADJ\n"
"Identify and pop simple bubbles.\n"
"\n"
"  -b, --bubble-length=N pop bubbles shorter than N bp\n"
"  -k, --kmer=K          pop bubbles shorter than 3*K bp\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static int k;
	static unsigned maxLength;
}

static const char shortopts[] = "b:k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bubble-length", required_argument, NULL, 'b' },
	{ "kmer",    required_argument, NULL, 'k' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

struct Contig {
	string id;
	unsigned length;
	unsigned coverage;
	Contig(const string& id, unsigned length, unsigned coverage)
		: id(id), length(length), coverage(coverage) { }
};

static vector<Contig> g_contigs;

typedef set<ContigNode> ContigNodes;
typedef map<ContigNode, ContigNodes> ContigGraph;
static ContigGraph g_graph;

inline unsigned ContigNode::outDegree() const
{
	return g_graph[*this].size();
}

inline unsigned ContigNode::inDegree() const
{
	return g_graph[~*this].size();
}

inline unsigned ContigNode::length() const
{
	return g_contigs[id()].length;
}

/** Return whether contig a has higher coverage than contig b. */
static bool compareCoverage(const ContigNode& a, const ContigNode& b)
{
	return g_contigs[a.id()].coverage > g_contigs[b.id()].coverage;
}

/** Popped branches. */
static set<unsigned> g_popped;

static void popBubble(const ContigNode& head,
		const ContigNodes& branches,
		const ContigNode& tail)
{
	assert(!branches.empty());
	assert(head.outDegree() == tail.inDegree());
	vector<ContigNode> sorted(branches.begin(), branches.end());
	sort(sorted.begin(), sorted.end(), compareCoverage);
	if (opt::verbose > 2) {
		cerr << '"' << head << "\" -> {";
		copy(sorted.begin(), sorted.end(),
				affix_ostream_iterator<ContigNode>(cerr,
					" \"", "\""));
		cerr << " } -> \"" << tail << "\";\n";
	}
	transform(sorted.begin() + 1, sorted.end(),
			inserter(g_popped, g_popped.begin()),
			mem_fun_ref(&ContigNode::id));
}

static struct {
	unsigned bubbles;
	unsigned popped;
	unsigned tooLong;
} g_count;

static void consider(const ContigNode& head,
		const ContigNodes& branches)
{
	assert(branches.size() > 1);
	ContigNodes tails;
	for (ContigNodes::const_iterator it = branches.begin();
			it != branches.end(); ++it) {
		const ContigNode& branch = *it;
		if (branch.outDegree() != 1 || branch.inDegree() != 1) {
			// This branch is not simple.
			return;
		}
		const ContigNode& tail = *g_graph[branch].begin();
		tails.insert(tail);
	}
	if (tails.size() != 1) {
		// The branches do not merge back to the same node.
		return;
	}

	const ContigNode& tail = *tails.begin();
	if (tail.inDegree() != branches.size()) {
		// This branch is not simple.
		return;
	}

	g_count.bubbles++;
	set<unsigned> lengths;
	transform(branches.begin(), branches.end(),
			inserter(lengths, lengths.begin()),
			mem_fun_ref(&ContigNode::length));
	unsigned minLength = *lengths.begin();
	unsigned maxLength = *lengths.rbegin();
	if (opt::verbose > 1)
		cerr << minLength << '\t' << maxLength << '\n';
	if (maxLength >= opt::maxLength) {
		// This branch is too long.
		g_count.tooLong++;
		return;
	}

	g_count.popped++;
	popBubble(head, branches, tail);
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

static void readContigGraph(ContigGraph& graph,
		const std::string& path)
{
	ifstream fin(path.c_str());
	istream& in = path.empty() || path == "-" ? cin : fin;
	if (&in == &fin)
		assert_open(fin, path);
	assert(in.good());

	string id;
	while (in >> id) {
		in.ignore(numeric_limits<streamsize>::max(), ';');
		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			ContigNode head(id, dir);
			string s;
			getline(in, s, dir == SENSE ? ';' : '\n');
			assert(in.good());
			istringstream ss(s);
			for (ContigNode tail; ss >> tail;) {
				if (dir == ANTISENSE)
					tail = ~tail;
				graph[head].insert(tail);
			}
			assert(ss.eof());
		}
	}

	assert(in.eof());
}

static void readContigs(vector<Contig>& contigs, const string& path)
{
	ifstream in(path.c_str());
	assert_open(in, path);
	assert(in.good());

	string id;
	unsigned length, coverage;
	while (in >> id >> length >> coverage) {
		in.ignore(numeric_limits<streamsize>::max(), '\n');
		unsigned serial = g_contigIDs.serial(id);
		assert(contigs.size() == serial);
		(void)serial;
		contigs.push_back(Contig(id, length, coverage));
	}
	assert(in.eof());
	assert(!contigs.empty());
	g_contigIDs.lock();
}

static const string& idToString(unsigned id)
{
	return g_contigIDs.key(id);
}

int main(int argc, char *const argv[])
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'b': arg >> opt::maxLength; break;
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

	if (opt::maxLength <= 0 && opt::k > 0)
		opt::maxLength = 3 * opt::k;

	if (opt::maxLength <= 0) {
		cerr << PROGRAM ": " << "missing -l,--length option\n";
		die = true;
	}

	if (argc - optind < 0) {
		cerr << PROGRAM ": missing arguments\n";
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

	string adjPath(optind < argc ? argv[optind++] : "-");
	readContigs(g_contigs, adjPath);
	readContigGraph(g_graph, adjPath);

	for (ContigGraph::const_iterator it = g_graph.begin();
			it != g_graph.end(); ++it) {
		const ContigGraph::value_type& node = *it;
		if (node.second.size() > 1)
			consider(node.first, node.second);
	}

	transform(g_popped.begin(), g_popped.end(),
			ostream_iterator<string>(cout, "\n"),
			idToString);

	if (opt::verbose > 0)
		cerr << "Bubbles: " << g_count.bubbles/2
			<< " Popped: " << g_popped.size()
			<< " Too long: " << g_count.tooLong/2
			<< '\n';
	return 0;
}
