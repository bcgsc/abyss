/**
 * Identify and pop simple bubbles.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "config.h"
#include "Common/Options.h"
#include "ConstString.h"
#include "ContigGraph.h"
#include "ContigGraphAlgorithms.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "DirectedGraph.h"
#include "FastaReader.h"
#include "GraphIO.h"
#include "Iterator.h"
#include "needleman_wunsch.h"
#include "Sequence.h"
#include "Uncompress.h"
#include <algorithm>
#include <cerrno>
#include <climits> // for UINT_MAX
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits> // for numeric_limits
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "PopBubbles"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... FASTA ADJ\n"
"Identify and pop simple bubbles.\n"
"  FASTA  contigs in FASTA format\n"
"  ADJ    contig adjacency graph\n"
"\n"
"  -k, --kmer=N          k-mer size\n"
"  -b, --bubble-length=N pop bubbles shorter than N bp\n"
"                        default is unlimited\n"
"  -p, --identity=REAL   minimum identity, default: 0.9\n"
"  -g, --graph=FILE      write the contig adjacency graph to FILE\n"
"      --dot             output bubbles in dot format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties

	/** Pop bubbles shorter than this threshold. */
	static unsigned maxLength = UINT_MAX;

	/** Minimum identity. */
	static float identity = 0.9;

	/** Write the contig adjacency graph to this file. */
	static string graphPath;

	/** Output bubbles in dot format. */
	static int dot;

	int format; // used by ContigProperties
}

static const char shortopts[] = "b:g:k:p:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bubble-length", required_argument, NULL, 'b' },
	{ "dot",           no_argument,       &opt::dot, 1, },
	{ "graph",         required_argument, NULL, 'g' },
	{ "kmer",          required_argument, NULL, 'k' },
	{ "identity",      required_argument, NULL, 'p' },
	{ "verbose",       no_argument,       NULL, 'v' },
	{ "help",          no_argument,       NULL, OPT_HELP },
	{ "version",       no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Popped branches. */
static vector<ContigID> g_popped;

/** Contig adjacency graph. */
typedef ContigGraph<DirectedGraph<ContigProperties, Distance> > Graph;
typedef Graph::vertex_descriptor vertex_descriptor;
typedef Graph::vertex_iterator vertex_iterator;
typedef Graph::adjacency_iterator adjacency_iterator;

struct CompareCoverage {
	const Graph& g;
	CompareCoverage(const Graph& g) : g(g) { }
	bool operator()(vertex_descriptor u, vertex_descriptor v)
	{
		return g[u].coverage > g[v].coverage;
	}
};

/** Pop the bubble between vertices v and tail. */
static void popBubble(Graph& g,
		vertex_descriptor v, vertex_descriptor tail)
{
	unsigned nbranches = g.out_degree(v);
	assert(nbranches > 1);
	assert(nbranches == g.in_degree(tail));
	vector<vertex_descriptor> sorted(nbranches);
	pair<adjacency_iterator, adjacency_iterator>
		adj = g.adjacent_vertices(v);
	copy(adj.first, adj.second, sorted.begin());
	sort(sorted.begin(), sorted.end(), CompareCoverage(g));
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
	unsigned tooMany;
	unsigned dissimilar;
} g_count;

/** Contig sequences. */
typedef vector<const_string> Contigs;
static Contigs g_contigs;

/** Return the sequence of vertex u. */
static string getSequence(vertex_descriptor u)
{
	assert(!u.ambiguous());
	assert(u.id() < g_contigs.size());
	string seq(g_contigs[u.id()]);
	return u.sense() ? reverseComplement(seq) : seq;
}

/** Return the length of vertex v. */
static unsigned getLength(const Graph* g, vertex_descriptor v)
{
	return (*g)[v].length;
}

/** Align the sequences of [first,last).
 * @return the identity of the global alignment
 */
template <typename It>
static float getAlignmentIdentity(It first, It last)
{
	assert(distance(first, last) == 2);
	(void)last;
	string seqa = getSequence(*first);
	++first;
	string seqb = getSequence(*first);

	NWAlignment alignment;
	unsigned matches = alignGlobal(seqa, seqb, alignment);
	if (opt::verbose > 2)
		cerr << alignment;
	return (float)matches / alignment.size();
}

/** Consider popping the bubble originating at the vertex v. */
static void considerPopping(Graph* pg, vertex_descriptor v)
{
	Graph& g = *pg;
	unsigned nbranches = g.out_degree(v);
	if (nbranches < 2)
		return;
	vertex_descriptor v1 = *g.adjacent_vertices(v).first;
	if (g.out_degree(v1) != 1) {
		// This branch is not simple.
		return;
	}
	vertex_descriptor tail = *g.adjacent_vertices(v1).first;
	if (v == ~tail // Palindrome
			|| g.in_degree(tail) != nbranches) {
		// This branch is not simple.
		return;
	}

	// Check that every branch is simple and ends at the same node.
	pair<adjacency_iterator, adjacency_iterator>
		adj = g.adjacent_vertices(v);
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

	if (opt::verbose > 2) {
		cerr << "\n* " << v << " -> ";
		copy(adj.first, adj.second,
				ostream_iterator<ContigNode>(cerr, " "));
		cerr << "-> " << tail << '\n';
	}

	g_count.bubbles++;
	const unsigned MAX_BRANCHES = opt::identity > 0 ? 2 : UINT_MAX;
	if (nbranches > MAX_BRANCHES) {
		// Too many branches.
		g_count.tooMany++;
		if (opt::verbose > 1)
			cerr << nbranches << " paths (too many)\n";
		return;
	}

	vector<unsigned> lengths(nbranches);
	transform(adj.first, adj.second, lengths.begin(),
			bind1st(ptr_fun(getLength), &g));
	unsigned minLength = *min_element(lengths.begin(), lengths.end());
	unsigned maxLength = *max_element(lengths.begin(), lengths.end());
	if (maxLength >= opt::maxLength) {
		// This branch is too long.
		g_count.tooLong++;
		if (opt::verbose > 1)
			cerr << minLength << '\t' << maxLength
				<< "\t0\t(too long)\n";
		return;
	}

	float identity = opt::identity == 0 ? 0
		: getAlignmentIdentity(adj.first, adj.second);
	bool dissimilar = identity < opt::identity;
	if (opt::verbose > 1)
		cerr << minLength << '\t' << maxLength << '\t' << identity
			<< (dissimilar ? "\t(dissimilar)" : "") << '\n';
	if (dissimilar) {
		// Insufficient identity.
		g_count.dissimilar++;
		return;
	}

	g_count.popped++;
	popBubble(g, v, tail);
}

/** Remove the specified contig from the adjacency graph. */
static void removeContig(Graph* g, ContigID id)
{
	ContigNode v(id, false);
	g->clear_vertex(v);
	g->remove_vertex(v);
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

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
			case 'b': arg >> opt::maxLength; break;
			case 'g': arg >> opt::graphPath; break;
			case 'k': arg >> opt::k; break;
			case 'p': arg >> opt::identity; break;
			case 'v': opt::verbose++; break;
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

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (argc - optind > 2) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	const char* contigsPath(argv[optind++]);
	string adjPath(argv[optind++]);

	// Read the contigs.
	Contigs& contigs = g_contigs;
	if (opt::identity > 0) {
		FastaReader in(contigsPath, FastaReader::NO_FOLD_CASE);
		for (FastaRecord rec; in >> rec;) {
			ContigID id(rec.id);
			assert(contigs.size() == id);
			contigs.push_back(rec.seq);
		}
		assert(in.eof());
		assert(!contigs.empty());
		opt::colourSpace = isdigit(contigs.front()[0]);
		ContigID::lock();
	}

	// Read the contig adjacency graph.
	ifstream fin(adjPath.c_str());
	assert_open(fin, adjPath);
	Graph g;
	fin >> g;
	assert(fin.eof());

	if (opt::dot)
		cout << "digraph bubbles {\n";
	pair<vertex_iterator, vertex_iterator> vit = g.vertices();
	for_each(vit.first, vit.second,
			bind1st(ptr_fun(considerPopping), &g));

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
			<< " Too many: " << g_count.tooMany/2
			<< " Dissimilar: " << g_count.dissimilar/2
			<< '\n';

	if (!opt::graphPath.empty()) {
		// Remove the popped contigs from the adjacency graph.
		for_each(g_popped.begin(), g_popped.end(),
				bind1st(ptr_fun(removeContig), &g));

		// Assemble unambiguous paths.
		typedef vector<ContigPath> ContigPaths;
		ContigPaths paths;
		assemble(g, back_inserter(paths));
		for (ContigPaths::const_iterator it = paths.begin();
				it != paths.end(); ++it)
			cout << ContigID::create() << '\t' << *it << '\n';
		paths.clear();

		// Output the updated adjacency graph.
		ofstream fout(opt::graphPath.c_str());
		assert(fout.good());
		write_graph(fout, g, PROGRAM, commandLine);
		assert(fout.good());
	}

	return 0;
}
