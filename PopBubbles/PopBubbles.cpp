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
#include "Sequence.h"
#include "Uncompress.h"
#include "alignGlobal.h"
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
#include <utility>
#include <vector>
#if _OPENMP
# include <omp.h>
#endif

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
"                        default is 10000\n"
"  -p, --identity=REAL   minimum identity, default: 0.9\n"
"      --scaffold        scaffold over bubbles that have\n"
"                        insufficient identity\n"
"      --no-scaffold     disable scaffolding [default]\n"
"  -g, --graph=FILE      write the contig adjacency graph to FILE\n"
"      --dot             output bubbles in dot format\n"
"  -j, --threads=N       use N parallel threads [1]\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties

	/** Pop bubbles shorter than this threshold. */
	static unsigned maxLength = 10000;

	/** Minimum identity. */
	static float identity = 0.9;

	/** Scaffold over bubbles that have insufficient identity. */
	static int scaffold;

	/** Write the contig adjacency graph to this file. */
	static string graphPath;

	/** Output bubbles in dot format. */
	static int dot;

	int format; // used by ContigProperties

	/** Number of threads. */
	static int threads = 1;
}

static const char shortopts[] = "b:g:j:k:p:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "bubble-length", required_argument, NULL, 'b' },
	{ "dot",           no_argument,       &opt::dot, 1, },
	{ "graph",         required_argument, NULL, 'g' },
	{ "kmer",          required_argument, NULL, 'k' },
	{ "identity",      required_argument, NULL, 'p' },
	{ "scaffold",      no_argument,       &opt::scaffold, 1},
	{ "no-scaffold",   no_argument,       &opt::scaffold, 0},
	{ "threads",       required_argument, NULL, 'j' },
	{ "verbose",       no_argument,       NULL, 'v' },
	{ "help",          no_argument,       NULL, OPT_HELP },
	{ "version",       no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Popped branches. */
static vector<ContigID> g_popped;

/** Bubbles that were not popped. */
static vector<pair<ContigNode, ContigNode> > g_bubbles;

/** Contig adjacency graph. */
typedef ContigGraph<DirectedGraph<ContigProperties, Distance> > Graph;
typedef Graph::vertex_descriptor vertex_descriptor;
typedef Graph::vertex_iterator vertex_iterator;
typedef Graph::adjacency_iterator adjacency_iterator;

/** Return the distance from vertex u to v. */
static int getDistance(const Graph& g,
		vertex_descriptor u, vertex_descriptor v)
{
	typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
	pair<edge_descriptor, bool> e = edge(u, v, g);
	assert(e.second);
	return g[e.first].distance;
}

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
	if (opt::dot)
#pragma omp critical(cout)
	{
		cout << '"' << v << "\" -> {";
		copy(sorted.begin(), sorted.end(),
				affix_ostream_iterator<ContigNode>(cout,
					" \"", "\""));
		cout << " } -> \"" << tail << "\"\n";
	}
#pragma omp critical(g_popped)
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
#pragma omp critical(cerr)
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

	if (opt::verbose > 2)
#pragma omp critical(cerr)
	{
		cerr << "\n* " << v << " -> ";
		copy(adj.first, adj.second,
				ostream_iterator<ContigNode>(cerr, " "));
		cerr << "-> " << tail << '\n';
	}

#pragma omp atomic
	g_count.bubbles++;
	const unsigned MAX_BRANCHES = opt::identity > 0 ? 2 : UINT_MAX;
	if (nbranches > MAX_BRANCHES) {
		// Too many branches.
#pragma omp atomic
		g_count.tooMany++;
		if (opt::verbose > 1)
#pragma omp critical(cerr)
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
#pragma omp atomic
		g_count.tooLong++;
		if (opt::verbose > 1)
#pragma omp critical(cerr)
			cerr << minLength << '\t' << maxLength
				<< "\t0\t(too long)\n";
		return;
	}

	float identity = opt::identity == 0 ? 0
		: getAlignmentIdentity(adj.first, adj.second);
	bool dissimilar = identity < opt::identity;
	if (opt::verbose > 1)
#pragma omp critical(cerr)
		cerr << minLength << '\t' << maxLength << '\t' << identity
			<< (dissimilar ? "\t(dissimilar)" : "") << '\n';
	if (dissimilar) {
		// Insufficient identity.
#pragma omp atomic
		g_count.dissimilar++;
		if (opt::scaffold) {
#pragma omp critical(g_bubbles)
			g_bubbles.push_back(make_pair(v, tail));
		}
		return;
	}

#pragma omp atomic
	g_count.popped++;
	popBubble(g, v, tail);
}

/** Add distances to a path. */
static ContigPath addDistance(const Graph& g, const ContigPath& path)
{
	ContigPath out;
	out.reserve(path.size());
	ContigNode u = path.front();
	out.push_back(u);
	for (ContigPath::const_iterator it = path.begin() + 1;
			it != path.end(); ++it) {
		ContigNode v = *it;
		int distance = getDistance(g, u, v);
		if (distance >= 0) {
			int numN = distance + opt::k - 1; // by convention
			assert(numN >= 0);
			numN = max(numN, 1);
			out.push_back(ContigNode(numN, 'N'));
		}
		out.push_back(v);
		u = v;
	}
	return out;
}

/** Scaffold over the bubble between vertices u and w. */
static void scaffoldBubble(Graph* pg,
		pair<vertex_descriptor, vertex_descriptor> uw)
{
	typedef graph_traits<Graph>::adjacency_iterator Ait;
	typedef graph_traits<Graph>::vertex_descriptor V;
	Graph& g = *pg;
	V u = uw.first, w = uw.second;
	if (edge(u, w, g).second) {
		// Already scaffolded.
		return;
	}

	pair<Ait, Ait> vrange = g.adjacent_vertices(u);
	g_popped.insert(g_popped.end(), vrange.first, vrange.second);
	int maxDistance = INT_MIN;
	for (Ait vit = vrange.first; vit != vrange.second; ++vit) {
		V v = *vit;
		int distance
			= getDistance(g, u, v)
			+ g[v].length
			+ getDistance(g, v, w);
		maxDistance = max(maxDistance, distance);
	}
	assert(maxDistance != INT_MIN);
	add_edge(u, w, max(maxDistance, 1), g);
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
			case 'j': arg >> opt::threads; break;
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
#if _OPENMP
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#pragma omp parallel
#pragma omp single
	for (vertex_iterator it = vit.first; it != vit.second; ++it)
#pragma omp task firstprivate(it)
		considerPopping(&g, *it);
#else
	for_each(vit.first, vit.second,
			bind1st(ptr_fun(considerPopping), &g));
#endif

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
			<< " Popped: " << g_count.popped/2
			<< " Too long: " << g_count.tooLong/2
			<< " Too many: " << g_count.tooMany/2
			<< " Dissimilar: " << g_count.dissimilar/2
			<< '\n';

	if (!opt::graphPath.empty()) {
		// Remove the popped contigs from the adjacency graph.
		for_each(g_popped.begin(), g_popped.end(),
				bind1st(ptr_fun(removeContig), &g));

		// Scaffold over the remaining bubbles.
		g_popped.clear();
		for_each(g_bubbles.begin(), g_bubbles.end(),
				bind1st(ptr_fun(scaffoldBubble), &g));
		sort(g_popped.begin(), g_popped.end());
		assert(unique(g_popped.begin(), g_popped.end())
				== g_popped.end());
		copy(g_popped.begin(), g_popped.end(),
				ostream_iterator<ContigID>(cout, "\n"));
		for_each(g_popped.begin(), g_popped.end(),
				bind1st(ptr_fun(removeContig), &g));

		// Assemble unambiguous paths.
		typedef vector<ContigPath> ContigPaths;
		ContigPaths paths;
		if (opt::scaffold) {
			Graph gorig = g;
			assemble(g, back_inserter(paths));
			for (ContigPaths::const_iterator it = paths.begin();
					it != paths.end(); ++it)
				cout << ContigID::create() << '\t'
					<< addDistance(gorig, *it) << '\n';
		} else {
			assemble(g, back_inserter(paths));
			for (ContigPaths::const_iterator it = paths.begin();
					it != paths.end(); ++it)
				cout << ContigID::create() << '\t' << *it << '\n';
		}
		paths.clear();

		// Output the updated adjacency graph.
		ofstream fout(opt::graphPath.c_str());
		assert(fout.good());
		write_graph(fout, g, PROGRAM, commandLine);
		assert(fout.good());
	}

	return 0;
}
