/**
 * Remove short contigs that do not contribute any relevant
 * information to the assembly.
 * Written by Tony Raymond <traymond@bcgsc.ca>
 */

#include "Common/Options.h"
#include "ContigID.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "FastaReader.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include <algorithm>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace rel_ops;
using boost::tie;

#define PROGRAM "abyss-filtergraph"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Tony Raymond.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... ADJ\n"
"Remove short contigs that do not contribute any relevant\n"
"information to the assembly.\n"
"  ADJ    contig adjacency graph\n"
"\n"
"  -k, --kmer=N          k-mer size\n"
"  -T, --island=N        remove islands shorter than N [0]\n"
"  -t, --tip=N           remove tips shorter than N [0]\n"
"  -l, --length=N        remove contigs shorter than N [0]\n"
"      --shim            remove filler contigs that only contribute\n"
"                        to adjacency\n"
"      --no-shim         disable filler contigs removal [default]\n"
"  -m, --min-overlap=N   require a minimum overlap of N bases [10]\n"
"      --assemble        assemble unambiguous paths\n"
"      --no-assemble     disable assembling of paths [default]\n"
"  -g, --graph=FILE      write the contig adjacency graph to FILE\n"
"  -i, --ignore=FILE     ignore contigs seen in FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties

	/** Remove island contigs less than this length. */
	static unsigned minIslandLen = 0;

	/** Remove tips less than this length. */
	static unsigned minTipLen = 0;

	/** Remove all contigs less than this length. */
	static unsigned minLen = 0;

	/** Remove short contigs that don't contribute any sequence. */
	static int shim = 0;

	/** Assemble unambiguous paths. */
	static int assemble = 0;

	/** Write the contig adjacency graph to this file. */
	static string graphPath;

	/** Contigs to ignore. */
	static string ignorePath;

	/** The minimum overlap allowed between two contigs. */
	static int minOverlap = 10;

	int format; // used by ContigProperties
}

static const char shortopts[] = "g:i:k:l:m:t:T:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "graph",         required_argument, NULL, 'g' },
	{ "ignore",        required_argument, NULL, 'i' },
	{ "kmer",          required_argument, NULL, 'k' },
	{ "island",        required_argument, NULL, 'T' },
	{ "tip",           required_argument, NULL, 't' },
	{ "length",        required_argument, NULL, 'l' },
	{ "shim",          no_argument,       &opt::shim, 1 },
	{ "no-shim",       no_argument,       &opt::shim, 0 },
	{ "assemble",      no_argument,       &opt::assemble, 1 },
	{ "no-assemble",   no_argument,       &opt::assemble, 0 },
	{ "min-overlap",   required_argument, NULL, 'm' },
	{ "verbose",       no_argument,       NULL, 'v' },
	{ "help",          no_argument,       NULL, OPT_HELP },
	{ "version",       no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static vector<ContigID> g_removed;

/** Contig adjacency graph. */
typedef ContigGraph<DirectedGraph<ContigProperties, Distance> > Graph;
typedef Graph::vertex_descriptor vertex_descriptor;

/** Data for verbose output. */
static struct {
	unsigned removed;
	unsigned tails;
	unsigned too_long;
	unsigned too_complex;
	unsigned self_adj;
	unsigned checked;
	unsigned parallel_edge;
} g_count;

/** Returns if the contig can be removed from the graph. */
static bool removable(const Graph* pg, vertex_descriptor v)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::out_edge_iterator OEit;
	typedef GTraits::in_edge_iterator IEit;

	const Graph& g = *pg;

	// Check if previously removed
	if (get(vertex_removed, g, v)) {
		g_count.removed++;
		return false;
	}

	// Check for tails
	if (out_degree(v, g) == 0 || in_degree(v, g) == 0) {
		g_count.tails++;
		return false;
	}

	// Check that the result will be less complex that the original
	if (!(out_degree(v, g) == 1 && in_degree(v, g) > 1)
			&& !(out_degree(v, g) > 1 && in_degree(v, g) == 1)) {
		g_count.too_complex++;
		return false;
	}

	// Check if self adjacent
	OEit oei0, oei1;
	tie(oei0, oei1) = out_edges(v, g);
	for (OEit vw = oei0; vw != oei1; ++vw) {
		if (v == target(*vw, g) || ~v == target(*vw, g)) {
			g_count.self_adj++;
			return false;
		}
	}

	// Check that removing the contig will result in adjacent contigs
	// overlapping by at least opt::minOverlap.
	IEit iei0, iei1;
	tie(iei0, iei1) = in_edges(v, g);
	IEit maxuv = iei0;
	for (IEit uv = iei0; uv != iei1; ++uv)
		if (g[*maxuv].distance < g[*uv].distance)
			maxuv = uv;
	OEit maxvw = oei0;
	for (OEit vw = oei0; vw != oei1; ++vw)
		if (g[*maxvw].distance < g[*vw].distance)
			maxvw = vw;

	if (g[*maxuv].distance + (int)g[v].length + g[*maxvw].distance >
			-opt::minOverlap) {
		g_count.too_long++;
		return false;
	}
	return true;
}

/** Data to store information of an edge. */
struct EdgeInfo {
	vertex_descriptor u;
	vertex_descriptor w;
	Distance ep;

	EdgeInfo(vertex_descriptor u, vertex_descriptor w, int ep)
		: u(u), w(w), ep(ep) {}
	EdgeInfo() : u(), w(), ep() {}
};

/** Returns a list of edges that may be added when the vertex v is
 * removed. */
static bool findNewEdges(const Graph& g, vertex_descriptor v,
		vector<EdgeInfo>& eds, vector<bool>& markedContigs)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_descriptor V;
	typedef GTraits::out_edge_iterator OEit;
	typedef GTraits::in_edge_iterator IEit;

	IEit iei0, iei1;
	tie(iei0, iei1) = in_edges(v, g);
	OEit oei0, oei1;
	tie(oei0, oei1) = out_edges(v, g);

	vector<V> marked;

	// if not marked and longest link LE contig length.
	// for every edge from u->v and v->w we must add an edge u->w
	for (IEit uv = iei0; uv != iei1; ++uv) {
		for (OEit vw = oei0; vw != oei1; ++vw) {
			int x = g[*uv].distance + (int)g[v].length +
				g[*vw].distance;
			assert(x <= 0);
			EdgeInfo ed(source(*uv, g), target(*vw, g), x);
			eds.push_back(ed);
			if (out_degree(v, g) > 1)
				marked.push_back(ed.u);
			if (in_degree(v, g) > 1)
				marked.push_back(ed.w);

			// Don't remove a vertex if the result is a parallel edge.
			if (edge(ed.u, ed.w, g).second) {
				g_count.parallel_edge++;
				return false;
			}
		}
	}
	for (vector<V>::const_iterator it = marked.begin(); it != marked.end();
			it++)
		markedContigs[get(vertex_index, g, *it)] = true;
	return true;
}

/** Adds all edges described in the vector eds. */
static void addNewEdges(Graph& g, const vector<EdgeInfo>& eds)
{
	for (vector<EdgeInfo>::const_iterator edsit = eds.begin();
			edsit != eds.end(); ++edsit) {
		assert(!edge(edsit->u, edsit->w, g).second);
		assert(edsit->ep.distance <= -opt::minOverlap);
		add_edge(edsit->u, edsit->w, edsit->ep, g);
	}
}

/** Remove the specified contig from the adjacency graph. */
static void removeContigs(Graph& g, vector<vertex_descriptor>& sc)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_descriptor V;
	typedef GTraits::out_edge_iterator OEit;
	typedef GTraits::in_edge_iterator IEit;

	vector<vertex_descriptor> out;
	out.reserve(sc.size());

	vector<bool> markedContigs(g.num_vertices());
	for (vector<vertex_descriptor>::iterator it = sc.begin();
			it != sc.end(); ++it) {
		V v = *it;
		if (opt::verbose > 0 && ++g_count.checked % 10000000 == 0)
			cerr << "Removed " << g_count.removed << "/"
				<< g_count.checked
				<< " vertices that have been checked.\n";

		if (markedContigs[get(vertex_index, g, v)]) {
			out.push_back(v);
			continue;
		}

		if (!removable(&g, v))
			continue;

		vector<EdgeInfo> eds;
		if (findNewEdges(g, v, eds, markedContigs))
			addNewEdges(g, eds);
		else
			continue;

		clear_vertex(v, g);
		remove_vertex(v, g);
		g_removed.push_back(v);
		g_count.removed++;
	}
	sc.swap(out);
}

/** Finds all potentially removable contigs in the graph. */
static void findShortContigs(const Graph& g,
		vector<vertex_descriptor>& sc)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_iterator Vit;
	Vit first, second;
	tie(first, second) = vertices(g);
	copy_if(first, second, back_inserter(sc),
			bind1st(ptr_fun(removable), &g));
}

/** Functor used for sorting contigs based on degree, then size,
 * and then ID. */
struct sortContigs {
	const Graph& g;

	sortContigs(const Graph& g) : g(g) {}

	template <typename V>
	bool operator() (V a, V b)
	{
		const ContigProperties& ap = g[a];
		const ContigProperties& bp = g[b];

		unsigned dega = out_degree(a, g) * in_degree(a, g);
		unsigned degb = out_degree(b, g) * in_degree(b, g);

		return dega != degb ? dega < degb
			: ap.length != bp.length ? ap.length < bp.length
			: a < b;
	}
};

struct ShorterThanX : unary_function<vertex_descriptor, bool> {
	const Graph& g;
	const vector<bool>& seen;
	size_t x;

	ShorterThanX(const Graph& g, const vector<bool>& seen, size_t x)
		: g(g), seen(seen), x(x) { }

	bool operator()(vertex_descriptor y) const
	{
		return g[y].length < x && !get(vertex_removed, g, y)
			&& !seen[ContigID(y)];
	}
};

static void removeShims(Graph& g)
{
	if (opt::verbose > 0)
		cerr << "Removing shims from graph...\n";
	vector<vertex_descriptor> shortContigs;
	findShortContigs(g, shortContigs);
	for (unsigned i = 0; !shortContigs.empty(); ++i) {
		if (opt::verbose > 0)
			cerr << "Starting pass " << i << ": There are "
				<< shortContigs.size() << " contigs to check.\n";
		sort(shortContigs.begin(), shortContigs.end(), sortContigs(g));
		removeContigs(g, shortContigs);
	}
	if (opt::verbose > 0) {
		cerr << "Shim removal stats:\n";
		cerr << "Removed: " << g_count.removed/2
			<< " Too Complex: " << g_count.too_complex/2
			<< " Tails: " << g_count.tails/2
			<< " Too Long: " << g_count.too_long/2
			<< " Self Adjacent: " << g_count.self_adj/2
			<< " Parallel Edges: " << g_count.parallel_edge/2 << '\n';
	}
}

static void removeShortContigs(Graph& g, const vector<bool>& seen)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_iterator Vit;
	typedef GTraits::vertex_descriptor V;
	Vit first, second;
	tie(first, second) = vertices(g);
	vector<V> sc;
	copy_if(first, second, back_inserter(sc),
			ShorterThanX(g, seen, opt::minLen));
	remove_vertex_if(g, sc.begin(), sc.end(), True<V>());
	copy(sc.begin(), sc.end(), back_inserter(g_removed));
	if (opt::verbose > 0)
		cerr << "Removed " << sc.size()/2 << " short contigs.\n";
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
			case 'l': arg >> opt::minLen; assert(arg.eof()); break;
			case 'm': arg >> opt::minOverlap; assert(arg.eof()); break;
			case 'g': arg >> opt::graphPath; assert(arg.eof()); break;
			case 'i': arg >> opt::ignorePath; assert(arg.eof()); break;
			case 'k': arg >> opt::k; assert(arg.eof()); break;
			case 'T': arg >> opt::minIslandLen; assert(arg.eof()); break;
			case 't': arg >> opt::minTipLen; assert(arg.eof()); break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (opt::minOverlap < 0) {
		cerr << PROGRAM ": "
			<< "--min-overlap must be a positive integer.\n";
		die = true;
	}

	if (opt::k <= 0) {
		cerr << PROGRAM ": " << "missing -k,--kmer option\n";
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

	Graph g;
	// Read the contig adjacency graph.
	{
		string adjPath(argv[optind++]);
		if (opt::verbose > 0)
			cerr << "Loading graph from file: " << adjPath
				<< '\n';
		ifstream fin(adjPath.c_str());
		assert_good(fin, adjPath);
		fin >> g;
		assert(fin.eof());
	}

	// Read the set of contigs to ignore.
	vector<bool> seen(num_vertices(g) / 2);
	if (!opt::ignorePath.empty()) {
		ifstream in(opt::ignorePath.c_str());
		assert_good(in, opt::ignorePath);
		markSeenInPath(in, seen);
	}

	if (opt::verbose > 0) {
		cerr << "Graph stats before:\n";
		printGraphStats(cerr, g);
	}

	// Remove islands.
	if (opt::minIslandLen > 0) {
		size_t s = g_removed.size();
		removeIslands_if(g, back_inserter(g_removed),
				ShorterThanX(g, seen, opt::minIslandLen));
		if (opt::verbose)
			cerr << "Removed " << g_removed.size() - s
				<< " islands.\n";
	}

	// Remove tips.
	if (opt::minTipLen > 0) {
		size_t s = g_removed.size();
		pruneTips_if(g, back_inserter(g_removed),
				ShorterThanX(g, seen, opt::minTipLen));
		if (opt::verbose)
			cerr << "Removed " << g_removed.size() - s
				<< " tips.\n";
	}

	// Remove short contigs.
	if (opt::minLen > 0)
		removeShortContigs(g, seen);

	// Remove shims.
	if (opt::shim)
		removeShims(g);

	if (opt::verbose > 0) {
		cerr << "Graph stats after:\n";
		printGraphStats(cerr, g);
	}

	sort(g_removed.begin(), g_removed.end());
	g_removed.erase(unique(g_removed.begin(), g_removed.end()),
			g_removed.end());
	copy(g_removed.begin(), g_removed.end(),
			ostream_iterator<ContigID>(cout, "\n"));

	// Assemble unambiguous paths.
	if (opt::assemble) {
		typedef vector<ContigPath> ContigPaths;
		ContigPaths paths;
		assemble(g, back_inserter(paths));
		for (ContigPaths::const_iterator it = paths.begin();
				it != paths.end(); ++it)
			cout << ContigID::create() << '\t' << *it << '\n';
	}

	// Output the updated adjacency graph.
	if (!opt::graphPath.empty()) {
		ofstream fout(opt::graphPath.c_str());
		assert_good(fout, opt::graphPath);
		write_graph(fout, g, PROGRAM, commandLine);
		assert_good(fout, opt::graphPath);
	}

	return 0;
}
