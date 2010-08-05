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
typedef Graph::edge_descriptor edge_descriptor;
typedef Graph::out_edge_iterator out_edge_iterator;

static void popBubble(vertex_iterator v, const ContigNode& tail)
{
	assert(v->out_degree() > 0);
	assert(v->out_degree() == g_graph.in_degree(tail));
	vector<ContigNode> sorted(v->out_degree());
	transform(v->begin(), v->end(), sorted.begin(), g_graph.target);
	sort(sorted.begin(), sorted.end(), compareCoverage);
	if (opt::dot) {
		cout << '"' << *v << "\" -> {";
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

/** Return the length of the target vertex of the specified edge. */
static unsigned targetLength(edge_descriptor e)
{
	return g_graph[g_graph.target(e)].length;
}

/** Consider popping the bubble originating at the vertex v. */
static void considerPopping(vertex_iterator v)
{
	const Graph& g = g_graph;
	assert(v->out_degree() > 1);
	if (g.out_degree(v->front().target()) != 1) {
		// This branch is not simple.
		return;
	}
	vertex_descriptor tail = g[v->front().target()].front().target();
	if (g.in_degree(tail) != v->out_degree()) {
		// This branch is not simple.
		return;
	}

	// Check that every branch is simple and ends at the same node.
	for (out_edge_iterator it = v->begin(); it != v->end(); ++it) {
		vertex_descriptor t = g.target(*it);
		if (g.out_degree(t) != 1 || g.in_degree(t) != 1) {
			// This branch is not simple.
			return;
		}
		if (g[t].front().target() != tail) {
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
	popBubble(v, tail);
}

/** Remove the specified contig from the adjacency graph. */
static void removeContig(ContigID id)
{
	ContigNode v(id, false);
	g_graph.clear_vertex(v);
	g_graph.remove_vertex(v);
}

typedef vector<vertex_descriptor> Path;

/** Return the sum of vp and the properties of vertex v. */
ContigProperties operator+(const ContigProperties& vp,
		vertex_descriptor v)
{
	assert(!v.ambiguous());
	return vp + g_graph[v]; 
}

/** Return the contig properties of the specified path. */
static ContigProperties calculateProperties(const Path& path)
{
	return accumulate(path.begin(), path.end(),
			ContigProperties(opt::k - 1, 0));
}

/** Merge the specified path and update the graph g. */
static void mergePath(Graph& g, const Path& path)
{
	ContigID id = ContigID::create();
	cout << id << '\t' << ContigPath(path) << '\n';
	vertex_descriptor v = g.add_vertex(calculateProperties(path));
	assert(ContigID(v) == id);
	g.copy_in_edges(path.front(), v);
	g.copy_out_edges(path.back(), v);
	for_each(path.begin(), path.end(), removeContig);
}

/** Assemble unambiguous paths. */
static void assemble(Graph& g)
{
	pair<vertex_iterator, vertex_iterator> vertices = g.vertices();
	for (vertex_iterator it = vertices.first;
			it != vertices.second; ++it) {
		if (g.contiguous_out(*it) && !g.contiguous_in(*it)) {
			Path path;
			assemble(g, *it, back_inserter(path));
			assert(path.size() >= 3);
			assert(path.front() != path.back());
			// Output only the canonical path.
			if (path.front() < path.back())
				mergePath(g, path);
		}
	}
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
	for (vertex_iterator it = g_graph.begin();
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
		assemble(g_graph);

		// Output the updated adjacency graph.
		ofstream fout(opt::graphPath.c_str());
		assert(fout.good());
		fout << g_graph;
		assert(fout.good());
	}

	return 0;
}
