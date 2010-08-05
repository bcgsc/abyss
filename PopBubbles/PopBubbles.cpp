/**
 * Identify and pop simple bubbles.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "config.h"
#include "Common/Options.h"
#include "AffixIterator.h"
#include "ContigGraph.h"
#include "ContigGraphAlgorithms.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "Sequence.h"
#include <algorithm>
#include <cerrno>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits> // for numeric_limits
#include <numeric> // for accumulate
#include <sstream>
#include <string>
#include <vector>

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
"  -g, --graph=FILE      write the contig adjacency graph to FILE\n"
"      --dot             output bubbles in dot format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	int k; // used by ContigLength
	static unsigned maxLength;

	/** Write the contig adjacency graph to this file. */
	static string graphPath;

	/** Output bubbles in dot format. */
	static int dot;
}

static const char shortopts[] = "b:g:k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bubble-length", required_argument, NULL, 'b' },
	{ "dot",           no_argument,       &opt::dot, 1, },
	{ "graph",         no_argument,       NULL, 'g' },
	{ "kmer",    required_argument, NULL, 'k' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Contig adjacency graph. */
typedef ContigGraph<ContigProperties> Graph;
static Graph g_graph;

/** Return whether contig a has higher coverage than contig b. */
static bool compareCoverage(const ContigNode& a, const ContigNode& b)
{
	return g_graph[a].coverage > g_graph[b].coverage;
}

/** Popped branches. */
static vector<ContigID> g_popped;

typedef Graph::vertex_descriptor vertex_descriptor;
typedef Graph::vertex_iterator vertex_iterator;
typedef Graph::adjacency_iterator adjacency_iterator;

/** Pop the bubble between vertices v and tail. */
static void popBubble(vertex_descriptor v, vertex_descriptor tail)
{
	unsigned nbranches = g_graph.out_degree(v);
	assert(nbranches > 1);
	assert(nbranches == g_graph.in_degree(tail));
	vector<vertex_descriptor> sorted(nbranches);
	pair<adjacency_iterator, adjacency_iterator>
		adj = g_graph.adjacent_vertices(v);
	copy(adj.first, adj.second, sorted.begin());
	sort(sorted.begin(), sorted.end(), compareCoverage);
	if (opt::dot) {
		cout << '"' << v << "\" -> {";
		copy(sorted.begin(), sorted.end(),
				affix_ostream_iterator<ContigNode>(cout,
					" \"", "\""));
		cout << " } -> \"" << tail << "\"\n";
	}
	transform(sorted.begin() + 1, sorted.end(),
			back_inserter(g_popped),
			mem_fun_ref(&ContigNode::operator ContigID));
}

static struct {
	unsigned bubbles;
	unsigned popped;
	unsigned tooLong;
} g_count;

/** Return the length of vertex v. */
static unsigned getLength(vertex_descriptor v)
{
	return g_graph[v].length;
}

/** Consider popping the bubble originating at the vertex v. */
static void considerPopping(vertex_descriptor v)
{
	const Graph& g = g_graph;
	unsigned nbranches = g.out_degree(v);
	if (nbranches < 2)
		return;
	vertex_descriptor v1 = *g.adjacent_vertices(v).first;
	if (g.out_degree(v1) != 1) {
		// This branch is not simple.
		return;
	}
	vertex_descriptor tail = *g.adjacent_vertices(v1).first;
	if (g.in_degree(tail) != nbranches) {
		// This branch is not simple.
		return;
	}

	// Check that every branch is simple and ends at the same node.
	pair<adjacency_iterator, adjacency_iterator>
		adj = g_graph.adjacent_vertices(v);
	for (adjacency_iterator it = adj.first; it != adj.second; ++it) {
		if (g.out_degree(*it) != 1 || g.in_degree(*it) != 1) {
			// This branch is not simple.
			return;
		}
		if (*g.adjacent_vertices(*it).first != tail) {
			// The branches do not merge back to the same node.
			return;
		}
	}

	g_count.bubbles++;
	vector<unsigned> lengths(nbranches);
	transform(adj.first, adj.second, lengths.begin(), getLength);
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
	popBubble(v, tail);
}

/** Remove the specified contig from the adjacency graph. */
static void removeContig(ContigID id)
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
			case 'g': arg >> opt::graphPath; break;
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
	pair<vertex_iterator, vertex_iterator> vit = g_graph.vertices();
	for_each(vit.first, vit.second, considerPopping);

	// Each bubble should be identified twice. Remove the duplicate.
	sort(g_popped.begin(), g_popped.end());
	g_popped.erase(unique(g_popped.begin(), g_popped.end()),
			g_popped.end());

	if (opt::dot)
		cout << "}\n";
	else
		copy(g_popped.begin(), g_popped.end(),
				ostream_iterator<ContigID>(cout, "\n"));

	if (opt::verbose > 0)
		cerr << "Bubbles: " << g_count.bubbles/2
			<< " Popped: " << g_popped.size()
			<< " Too long: " << g_count.tooLong/2
			<< '\n';

	if (!opt::graphPath.empty()) {
		// Remove the popped contigs from the adjacency graph.
		for_each(g_popped.begin(), g_popped.end(), removeContig);

		// Assemble unambiguous paths.
		assemble(g_graph, cout);

		// Output the updated adjacency graph.
		ofstream fout(opt::graphPath.c_str());
		assert(fout.good());
		fout << g_graph;
		assert(fout.good());
	}

	return 0;
}
