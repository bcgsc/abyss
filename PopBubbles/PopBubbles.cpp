/**
 * Identify and pop simple bubbles.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "config.h"
#include "Common/Options.h"
#include "AffixIterator.h"
#include "ContigGraph.h"
#include "ContigNode.h"
#include "Sequence.h"
#include <algorithm>
#include <cerrno>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits> // for numeric_limits
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
"      --dot             output bubbles in dot format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	int k; // used by ContigLength
	static unsigned maxLength;

	/** Output bubbles in dot format. */
	static int dot;
}

static const char shortopts[] = "b:k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bubble-length", required_argument, NULL, 'b' },
	{ "dot",           no_argument,       &opt::dot, 1, },
	{ "kmer",    required_argument, NULL, 'k' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Vertex attribute. */
struct Contig {
	string id;
	unsigned length;
	unsigned coverage;
	Contig(const string& id, unsigned length, unsigned coverage)
		: id(id), length(length), coverage(coverage) { }
};

/** Vertex attributes. */
static vector<Contig> g_contigs;

/** Contig adjacency graph. */
static ContigGraph g_graph;

/** Collection of edges. */
typedef ContigGraph::Edges Edges;

inline unsigned ContigNode::outDegree() const
{
	return g_graph[*this].out_degree();
}

/** Return whether contig a has higher coverage than contig b. */
static bool compareCoverage(const ContigNode& a, const ContigNode& b)
{
	return g_contigs[a.id()].coverage > g_contigs[b.id()].coverage;
}

/** Popped branches. */
static vector<unsigned> g_popped;

static void popBubble(const ContigNode& head, const Edges& branches,
		const ContigNode& tail)
{
	assert(!branches.empty());
	assert(g_graph.out_degree(head) == g_graph.in_degree(tail));
	vector<ContigNode> sorted(branches.size());
	transform(branches.begin(), branches.end(), sorted.begin(),
			mem_fun_ref(&Edges::value_type::target_key));
	sort(sorted.begin(), sorted.end(), compareCoverage);
	if (opt::dot) {
		cout << '"' << head << "\" -> {";
		copy(sorted.begin(), sorted.end(),
				affix_ostream_iterator<ContigNode>(cout,
					" \"", "\""));
		cout << " } -> \"" << tail << "\"\n";
	}
	transform(sorted.begin() + 1, sorted.end(),
			back_inserter(g_popped),
			mem_fun_ref(&ContigNode::id));
}

static struct {
	unsigned bubbles;
	unsigned popped;
	unsigned tooLong;
} g_count;

/** Return the length of the target of the specified edge. */
static unsigned targetLength(const Edges::value_type& e)
{
	return e.target_key().length();
}

static void consider(const ContigNode& head, const Edges& branches)
{
	assert(branches.size() > 1);
	set<ContigNode> tails;
	if (branches.front().target().out_degree() != 1) {
		// This branch is not simple.
		return;
	}
	const ContigNode& tail = branches.front().target()
		.out_edges().front().target_key();
	if (g_graph.in_degree(tail) != branches.size()) {
		// This branch is not simple.
		return;
	}

	// Check that every branch is simple and ends at the same node.
	for (Edges::const_iterator it = branches.begin();
			it != branches.end(); ++it) {
		if (it->target().out_degree() != 1
				|| g_graph.in_degree(it->target_key()) != 1) {
			// This branch is not simple.
			return;
		}
		if (it->target().out_edges().begin()->target_key() != tail) {
			// The branches do not merge back to the same node.
			return;
		}
	}

	g_count.bubbles++;
	set<unsigned> lengths;
	transform(branches.begin(), branches.end(),
			inserter(lengths, lengths.begin()),
			targetLength);
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

/** Read the contig attributes. */
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
		cerr << PROGRAM ": " << "missing -b,--bubble-length option\n";
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
	readContigGraph(g_graph, adjPath);
	g_contigIDs.lock();
	g_contigs.reserve(g_contigIDs.size());
	readContigs(g_contigs, adjPath);

	if (opt::dot)
		cout << "digraph bubbles {\n";
	for (ContigGraph::const_iterator it = g_graph.begin();
			it != g_graph.end(); ++it) {
		if (it->out_degree() > 1)
			consider(it->vertex(), it->out_edges());
	}

	// Each bubble should be identified twice. Remove the duplicate.
	sort(g_popped.begin(), g_popped.end());
	g_popped.erase(unique(g_popped.begin(), g_popped.end()),
			g_popped.end());

	if (opt::dot)
		cout << "}\n";
	else
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
