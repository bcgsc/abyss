/**
 * Remove short contigs that do not contribute any relevant
 * information to the assembly.
 * Written by Tony Raymond <traymond@bcgsc.ca>
 */

#include "Common/Options.h"
#include "ContigGraph.h"
#include "ContigGraphAlgorithms.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "DirectedGraph.h"
#include "FastaReader.h"
#include "GraphIO.h"
#include "GraphUtil.h"
#include "IOUtil.h"
#include "Uncompress.h"
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
"  -g, --graph=FILE      write the contig adjacency graph to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties

	/** Write the contig adjacency graph to this file. */
	static string graphPath;

	/** The maximum distance allowed for a new edge. */
	static int maxDistance = 0;

	int format; // used by ContigProperties
}

static const char shortopts[] = "d:g:k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "graph",         required_argument, NULL, 'g' },
	{ "kmer",          required_argument, NULL, 'k' },
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
	unsigned prev_removed;
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
		g_count.prev_removed++;
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
	// overlapping by at least opt::maxDistance.
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
			opt::maxDistance) {
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
		vector<EdgeInfo>& eds,
		vector<bool>& markedContigs)
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
		assert(edsit->ep.distance <= opt::maxDistance);
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

	vector<bool> markedContigs;
	for (vector<vertex_descriptor>::iterator it = sc.begin();
			it != sc.end(); ++it) {
		V v = *it;
		if (opt::verbose > 0 && ++g_count.checked % 10000000 == 0)
			cerr << "Removed " << g_count.removed << "/"
				<< g_count.checked
				<< " contigs that have been checked.\n";

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
			case 'd': arg >> opt::maxDistance; assert(arg.eof()); break;
			case 'g': arg >> opt::graphPath; assert(arg.eof()); break;
			case 'k': arg >> opt::k; assert(arg.eof()); break;
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

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (argc - optind > 1) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (opt::graphPath.empty()) {
		cerr << PROGRAM ": missing -g option\n";
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

	if (opt::verbose > 0) {
		cerr << "Graph stats before:\n";
		printGraphStats(cerr, g);
	}

	vector<vertex_descriptor> shortContigs;
	findShortContigs(g, shortContigs);
	for (unsigned i = 0; !shortContigs.empty(); ++i) {
		if (opt::verbose > 0)
			cerr << "Starting pass " << i << ": There are "
				<< shortContigs.size() << " contigs to check.\n";
		sort(shortContigs.begin(), shortContigs.end(), sortContigs(g));
		removeContigs(g, shortContigs);
		if (opt::verbose > 0)
			cerr << "Removed " << g_count.removed << " contigs.\n";
	}

	if (opt::verbose > 0) {
		cerr << "Graph stats after:\n";
		printGraphStats(cerr, g);
		cerr << "Removed: " << g_count.removed
			<< " Too Complex: " << g_count.too_complex
			<< " Tails: " << g_count.tails
			<< " Previously Removed: " << g_count.prev_removed
			<< " Too Long: " << g_count.too_long
			<< " Self Adjacent: " << g_count.self_adj
			<< " Parallel Edges: " << g_count.parallel_edge << '\n';
	}

	sort(g_removed.begin(), g_removed.end());
	copy(g_removed.begin(), g_removed.end(),
			ostream_iterator<ContigID>(cout, "\n"));

	// Output the updated adjacency graph.
	ofstream fout(opt::graphPath.c_str());
	assert_good(fout, opt::graphPath);
	write_graph(fout, g, PROGRAM, commandLine);
	assert_good(fout, opt::graphPath);

	return 0;
}
