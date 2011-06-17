#include "config.h"
#include "ContigGraph.h"
#include "ContigGraphAlgorithms.h"
#include "ContigNode.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "DirectedGraph.h"
#include "DotIO.h"
#include "Estimate.h"
#include "GraphAlgorithms.h"
#include "GraphUtil.h"
#include "IOUtil.h"
#include <cassert>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

using namespace std;

#define PROGRAM "abyss-scaffold"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [DIST]\n"
"Scaffold contigs using the distance estimate graph.\n"
"  DIST  estimates of the distance between contigs\n"
"\n"
"  -n, --npairs=N        minimum number of pairs [0]\n"
"  -s, --seed-length=N   minimum contig length [0]\n"
"  -k, --kmer=N          length of a k-mer\n"
"      --min-gap=N       minimum scaffold gap length to output\n"
"  -o, --out=FILE        write the paths to FILE\n"
"  -g, --graph=FILE      write the graph to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties
	int dot = true; // used by Estimate

	/** Minimum number of pairs. */
	static unsigned minNumPairs;

	/** Minimum contig length. */
	static unsigned minContigLength;

	/** Minimum scaffold gap length to output. */
	static int minGap = INT_MIN;

	/** Write the paths to this file. */
	static string out;

	/** Write the graph to this file. */
	static string graphPath;

	/** Verbose output. */
	static int verbose;
}

static const char shortopts[] = "g:k:n:o:s:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MIN_GAP };

static const struct option longopts[] = {
	{ "graph",       no_argument,       NULL, 'g' },
	{ "kmer",        required_argument, NULL, 'k' },
	{ "min-gap",     required_argument, NULL, OPT_MIN_GAP },
	{ "npairs",      required_argument, NULL, 'n' },
	{ "out",         required_argument, NULL, 'o' },
	{ "seed-length", required_argument, NULL, 's' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Contig length property. */
struct Length {
	unsigned length;

	Length& operator+=(const Length& o)
	{
		length += o.length;
		return *this;
	}

	template <typename T>
	Length& operator+=(const T& o)
	{
		assert((int)length + (int)o.distance > 0);
		length += o.distance;
		return *this;
	}

	friend ostream& operator<<(ostream& out, const Length& o)
	{
		return out << "l=" << o.length;
	}

	friend istream& operator>>(istream& in, Length& o)
	{
		return in >> expect("l =") >> o.length;
	}
};

/** A distance estimate graph. */
typedef DirectedGraph<Length, DistanceEst> DG;
typedef ContigGraph<DG> Graph;

/** Add missing complementary edges. */
static void addComplementaryEdges(DG& g)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::edge_descriptor E;
	typedef GTraits::edge_iterator Eit;
	typedef GTraits::vertex_descriptor V;

	std::pair<Eit, Eit> erange = edges(g);
	unsigned numAdded = 0;
	for (Eit eit = erange.first; eit != erange.second; ++eit) {
		E e = *eit;
		V u = source(e, g), v = target(e, g);
		if (!edge(~v, ~u, g).second) {
			add_edge(~v, ~u, g[e], g);
			numAdded++;
		}
	}
	if (opt::verbose > 0)
		cerr << "Added " << numAdded << " complementary edges.\n";
}

/** Return whether the specified edges has sufficient support. */
struct PoorSupport {
	PoorSupport(Graph& g) : m_g(g) { }
	bool operator()(graph_traits<Graph>::edge_descriptor e) const
	{
		return m_g[e].numPairs < opt::minNumPairs;
	}
	const Graph& m_g;
};

/** Remove short vertices and unsupported edges from the graph. */
static void filterGraph(Graph& g)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::edge_descriptor E;
	typedef GTraits::edge_iterator Eit;
	typedef GTraits::vertex_descriptor V;
	typedef GTraits::vertex_iterator Vit;

	// Remove short contigs.
	unsigned numRemovedV = 0;
	std::pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		if (g[u].length < opt::minContigLength)
			clear_vertex(u, g);
		if (out_degree(u, g) == 0 && in_degree(u, g) == 0) {
			remove_vertex(u, g);
			numRemovedV++;
		}
	}
	if (opt::verbose > 0)
		cerr << "Removed " << numRemovedV << " vertices.\n";

	// Remove poorly-supported edges.
	unsigned numBefore = num_edges(g);
	remove_edge_if(PoorSupport(g), static_cast<DG&>(g));
	unsigned numRemovedE = numBefore - num_edges(g);
	if (opt::verbose > 0)
		cerr << "Removed " << numRemovedE << " edges.\n";
}

/** Find edges in g0 that resolve forks in g.
 * For a pair of edges (u,v1) and (u,v2) in g, if exactly one of the
 * edges (v1,v2) or (v2,v1) exists in g0, add that edge to g.
 */
static void resolveForks(Graph& g, const Graph& g0)
{
	typedef graph_traits<Graph>::adjacency_iterator Vit;
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef graph_traits<Graph>::vertex_iterator Uit;
	typedef graph_traits<Graph>::vertex_descriptor V;

	unsigned numEdges = 0;
	pair<Uit, Uit> urange = vertices(g);
	for (Uit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		if (out_degree(u, g) < 2)
			continue;
		pair<Vit, Vit> vrange = adjacent_vertices(u, g);
		for (Vit vit1 = vrange.first; vit1 != vrange.second;) {
			V v1 = *vit1;
			++vit1;
			assert(v1 != u);
			for (Vit vit2 = vit1; vit2 != vrange.second; ++vit2) {
				V v2 = *vit2;
				assert(v2 != u);
				assert(v1 != v2);
				if (edge(v1, v2, g).second || edge(v2, v1, g).second)
					continue;
				pair<E, bool> e12 = edge(v1, v2, g0);
				pair<E, bool> e21 = edge(v2, v1, g0);
				if (e12.second && e21.second) {
					if (opt::verbose > 1)
						cerr << "cycle: " << v1 << ' ' << v2 << '\n';
				} else if (e12.second || e21.second) {
					E e = e12.second ? e12.first : e21.first;
					V v = source(e, g0), w = target(e, g0);
					add_edge(v, w, g0[e], g);
					numEdges++;
					if (opt::verbose > 1)
						cerr << v << " -> " << w
							<< " [" << g0[e] << "]\n";
				}
			}
		}
	}
	if (opt::verbose > 0)
		cerr << "Added " << numEdges
			<< " edges to ambiguous vertices.\n";
}

/** Remove tips.
 * For an edge (u,v), remove the vertex v if deg+(u) > 1
 * and deg-(v) = 1 and deg+(v) = 0.
 */
static void pruneTips(Graph& g)
{
	typedef graph_traits<Graph>::adjacency_iterator Vit;
	typedef graph_traits<Graph>::vertex_iterator Uit;
	typedef graph_traits<Graph>::vertex_descriptor V;

	/** Identify the tips. */
	vector<V> tips;
	pair<Uit, Uit> urange = vertices(g);
	for (Uit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		if (out_degree(u, g) < 2)
			continue;
		pair<Vit, Vit> vrange = adjacent_vertices(u, g);
		for (Vit vit = vrange.first; vit != vrange.second; ++vit) {
			V v = *vit;
			assert(v != u);
			if (in_degree(v, g) == 1 && out_degree(v, g) == 0)
				tips.push_back(v);
		}
	}

	/** Remove the tips. */
	for (vector<V>::const_iterator it = tips.begin();
			it != tips.end(); ++it) {
		V u = *it;
		clear_vertex(u, g);
		remove_vertex(u, g);
	}

	if (opt::verbose > 0) {
		cerr << "Removed " << tips.size() << " tips.\n";
		printGraphStats(cerr, g);
	}
}

/** Return whether the specified distance estimate is an exact
 * overlap.
 */
static bool isOverlap(const DistanceEst& d)
{
	if (d.stdDev == 0 && d.numPairs == 0) {
		assert(d.distance < 0);
		return true;
	} else
		return false;
}

/** Add distance estimates to a path. */
static ContigPath addDistEst(const Graph& g, const ContigPath& path)
{
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef edge_bundle_type<Graph>::type EP;

	ContigPath out;
	out.reserve(2 * path.size());
	ContigNode u = path.front();
	out.push_back(u);
	for (ContigPath::const_iterator it = path.begin() + 1;
			it != path.end(); ++it) {
		ContigNode v = *it;
		assert(!v.ambiguous());
		pair<E, bool> e = edge(u, v, g);
		assert(e.second);
		const EP& ep = g[e.first];
		if (!isOverlap(ep)) {
			int distance = max(ep.distance, (int)opt::minGap);
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

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'g': arg >> opt::graphPath; break;
			case 'n': arg >> opt::minNumPairs; break;
			case 'o': arg >> opt::out; break;
			case 's': arg >> opt::minContigLength; break;
			case 'v': opt::verbose++; break;
			case OPT_MIN_GAP: arg >> opt::minGap; break;
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

	if (argc - optind < 0) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	} else if (argc - optind > 1) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	string distFilePath(optind < argc ? argv[optind++] : "-");

	// Read the distance estimate graph.
	ifstream fin(distFilePath.c_str());
	istream& in = distFilePath == "-" ? cin : fin;
	assert_good(in, distFilePath);
	Graph g;
	read_dot<DG>(in, g);
	assert(in.eof());
	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	// Add any missing complementary edges.
	addComplementaryEdges(g);
	Graph gorig(g);

	// Filter the graph.
	filterGraph(g);
	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	// Resolve forks.
	resolveForks(g, gorig);

	// Prune tips.
	pruneTips(g);

	// Remove transitive edges.
	unsigned numTransitive = remove_transitive_edges(g);
	if (opt::verbose > 0) {
		cerr << "Removed " << numTransitive << " transitive edges.\n";
		printGraphStats(cerr, g);
	}

	// Assemble the paths.
	typedef vector<ContigPath> ContigPaths;
	ContigPaths paths;
	assemble(g, back_inserter(paths));
	sort(paths.begin(), paths.end());
	if (opt::verbose > 0) {
		unsigned n = 0;
		for (ContigPaths::const_iterator it = paths.begin();
				it != paths.end(); ++it)
			n += it->size();
		cerr << "Assembled " << n << " contigs in "
			<< paths.size() << " scaffolds.\n";
		printGraphStats(cerr, g);
	}

	// Output the paths.
	ofstream fout(opt::out.c_str());
	ostream& out = opt::out.empty() || opt::out == "-" ? cout : fout;
	assert_good(out, opt::out);
	for (vector<ContigPath>::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		out << ContigID::create() << '\t'
			<< addDistEst(gorig, *it) << '\n';
	}
	assert_good(out, opt::out);

	// Output the graph.
	if (!opt::graphPath.empty()) {
		ofstream out(opt::graphPath.c_str());
		assert_good(out, opt::graphPath);
		write_dot(out, g);
		assert_good(out, opt::graphPath);
	}

	return 0;
}
