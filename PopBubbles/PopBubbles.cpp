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

/** Contig properties. */
struct Contig {
	unsigned length;
	unsigned coverage;

	friend ostream& operator <<(ostream& out, const Contig& o)
	{
		return out << ' ' << o.length << ' ' << o.coverage;
	}

	friend istream& operator >>(istream& in, Contig& o)
	{
		return in >> o.length >> o.coverage;
	}
};

/** Contig adjacency graph. */
typedef ContigGraph<Contig> Graph;
static Graph g_graph;

/** Return whether contig a has higher coverage than contig b. */
static bool compareCoverage(const ContigNode& a, const ContigNode& b)
{
	return g_graph[a].coverage > g_graph[b].coverage;
}

/** Popped branches. */
static vector<unsigned> g_popped;

/** Return the target vertex of edge e. */
static ContigNode target(const Graph::Edge& e)
{
	return g_graph.target(e);
}

static void popBubble(Graph::vertex_iterator v,
		const ContigNode& tail)
{
	assert(v->out_degree() > 0);
	assert(v->out_degree() == g_graph.in_degree(tail));
	vector<ContigNode> sorted(v->out_degree());
	transform(v->begin(), v->end(), sorted.begin(), target);
	sort(sorted.begin(), sorted.end(), compareCoverage);
	if (opt::dot) {
		cout << '"' << g_graph.vertex(*v) << "\" -> {";
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

/** Return the length of the target vertex of the specified edge. */
static unsigned targetLength(const Graph::Edge& e)
{
	return g_graph[g_graph.target(e)].length;
}

/** Consider popping the bubble originating at the vertex v. */
static void considerPopping(Graph::vertex_iterator v)
{
	assert(v->out_degree() > 1);
	if (v->front().target().out_degree() != 1) {
		// This branch is not simple.
		return;
	}
	const Graph::Vertex& tail = v->front().target().front().target();
	if (g_graph.in_degree(tail) != v->out_degree()) {
		// This branch is not simple.
		return;
	}

	// Check that every branch is simple and ends at the same node.
	for (Graph::out_edge_iterator it = v->begin();
			it != v->end(); ++it) {
		if (it->target().out_degree() != 1
				|| g_graph.in_degree(it->target()) != 1) {
			// This branch is not simple.
			return;
		}
		if (it->target().front().target() != tail) {
			// The branches do not merge back to the same node.
			return;
		}
	}

	g_count.bubbles++;
	vector<unsigned> lengths(v->out_degree());
	transform(v->begin(), v->end(), lengths.begin(), targetLength);
	unsigned minLength = *min_element(lengths.begin(), lengths.end());
	unsigned maxLength = *max_element(lengths.begin(), lengths.end());
	if (opt::verbose > 1)
		cerr << minLength << '\t' << maxLength << '\n';
	if (maxLength >= opt::maxLength) {
		// This branch is too long.
		g_count.tooLong++;
		return;
	}

	g_count.popped++;
	popBubble(v, g_graph.vertex(tail));
}

/** Remove the specified contig from the adjacency graph. */
static void removeContig(unsigned id)
{
	ContigNode v(id, false);
	g_graph.clear_vertex(v);
	g_graph.remove_vertex(v);
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
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

	if (argc - optind < 1) {
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

	string adjPath(argv[optind++]);
	ifstream fin(adjPath.c_str());
	assert_open(fin, adjPath);
	fin >> g_graph;
	assert(fin.eof());

	if (opt::dot)
		cout << "digraph bubbles {\n";
	for (Graph::vertex_iterator it = g_graph.begin();
			it != g_graph.end(); ++it)
		if (it->out_degree() > 1)
			considerPopping(it);

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

	if (opt::verbose < 3)
		return 0;

	// Remove the popped contigs from the adjacency graph.
	for_each(g_popped.begin(), g_popped.end(), removeContig);

	// Output the updated adjacency graph.
	cerr << g_graph;
	assert(cerr.good());

	return 0;
}
