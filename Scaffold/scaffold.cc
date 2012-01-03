#include "config.h"
#include "ContigNode.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "Estimate.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphAlgorithms.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include "Graph/PopBubbles.h"
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
using namespace std::rel_ops;
using boost::edge_bundle_type;
using boost::tie;

#define PROGRAM "abyss-scaffold"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... FASTA|OVERLAP DIST...\n"
"Scaffold contigs using the distance estimate graph.\n"
"  FASTA    contigs in FASTA format\n"
"  OVERLAP  the contig overlap graph\n"
"  DIST     estimates of the distance between contigs\n"
"\n"
"  -n, --npairs=N        minimum number of pairs [0]\n"
"  -s, --seed-length=N   minimum contig length [0]\n"
"  -k, --kmer=N          length of a k-mer\n"
"      --min-gap=N       minimum scaffold gap length to output [50]\n"
"  -o, --out=FILE        write the paths to FILE\n"
"  -g, --graph=FILE      write the graph to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties

	/** Minimum number of pairs. */
	static unsigned minNumPairs;

	/** Minimum contig length. */
	static unsigned minContigLength;

	/** Minimum scaffold gap length to output. */
	static int minGap = 50;

	/** Write the paths to this file. */
	static string out;

	/** Write the graph to this file. */
	static string graphPath;

	/** Verbose output. */
	int verbose; // used by PopBubbles

 	/** Output format */
 	int format = DOT; // used by DistanceEst
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

	bool operator==(const Length& o) const
	{
		return length == o.length;
	}

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
		if (in >> std::ws && in.peek() == 'l')
			return in >> expect("l =") >> o.length;
		else
			return in >> o.length;
	}
};

static inline
void put(vertex_length_t, Length& vp, unsigned length)
{
	vp.length = length;
}

static inline
void put(vertex_coverage_t, Length&, unsigned)
{
}

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
						cerr << u << " -> " << v << " -> " << w
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
	typedef graph_traits<Graph>::vertex_descriptor V;

	/** Identify the tips. */
	vector<V> tips;
	pruneTips(g, back_inserter(tips));

	if (opt::verbose > 0) {
		cerr << "Removed " << tips.size() << " tips.\n";
		printGraphStats(cerr, g);
	}
}

/** Return true if the specified vertex is a repeat. */
static bool isRepeat(const Graph& g,
		graph_traits<Graph>::vertex_descriptor t1,
		graph_traits<Graph>::vertex_descriptor t2,
		graph_traits<Graph>::vertex_descriptor u,
		graph_traits<Graph>::vertex_descriptor v1,
		graph_traits<Graph>::vertex_descriptor v2)
{
	assert(in_degree(u, g) == 2 && out_degree(u, g) == 2);
	(void)u;
	return out_degree(t1, g) == 2 && out_degree(t2, g) == 2
		&& in_degree(v1, g) == 2 && in_degree(v2, g) == 2
		&& edge(t1, v1, g).second && edge(t2, v2, g).second
		&& !edge(t1, v2, g).second && !edge(t2, v1, g).second;
}

/** Remove repetitive vertices from this graph.
 * input: digraph g { t1->v1 t2->v2 t1->u t2->u u->v1 u->v2 }
 * operation: remove vertex u
 * output: digraph g { t1->v1 t2->v2 }
 */
static void removeRepeats(Graph& g)
{
	typedef graph_traits<Graph>::adjacency_iterator Ait;
	typedef graph_traits<Graph>::vertex_iterator Vit;
	typedef graph_traits<Graph>::vertex_descriptor V;

	vector<V> repeats;
	pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		if (in_degree(u, g) != 2 || out_degree(u, g) != 2)
			continue;

		Ait vit, vlast;
		tie(vit, vlast) = adjacent_vertices(u, g);
		V v1 = *vit;
		++vit;
		V v2 = *vit;
		++vit;
		assert(vit == vlast);
		assert(v1 != v2 && v1 != u && v2 != u);

		Ait tit, tlast;
		tie(tit, tlast) = adjacent_vertices(~u, g);
		V t1 = ~*tit;
		++tit;
		V t2 = ~*tit;
		++tit;
		assert(tit == tlast);
		assert(t1 != t2 && t1 != u && t2 != u);

		if (isRepeat(g, t1, t2, u, v1, v2))
			repeats.push_back(u);
		else if (isRepeat(g, t1, t2, u, v2, v1))
			repeats.push_back(u);
	}

	/** Remove the repetitive vertices. */
	for (vector<V>::const_iterator it = repeats.begin();
			it != repeats.end(); ++it) {
		V u = *it;
		clear_vertex(u, g);
		remove_vertex(u, g);
	}
	if (opt::verbose > 0) {
		cerr << "Removed "
			<< repeats.size() << " repetitive vertices.\n";
		printGraphStats(cerr, g);
	}
	if (opt::verbose > 1) {
		cerr << "Repeats: ";
		copy(repeats.begin(), repeats.end(),
				ostream_iterator<ContigNode>(cerr, " "));
		cerr << '\n';
	}
}

/** Remove weak edges from this graph.
 * input: digraph g { u1->v2 u1->v1 u2->v2 }
 *        (u1,v2).n < (u1,v1).n and (u1,v2).n < (u2,v2).n
 * operation: remove edge u1->v2
 * output: digraph g {u1->v1 u2->v2 }
 */
static void removeWeakEdges(Graph& g)
{
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef graph_traits<Graph>::edge_iterator Eit;
	typedef graph_traits<Graph>::in_edge_iterator Iit;
	typedef graph_traits<Graph>::out_edge_iterator Oit;
	typedef graph_traits<Graph>::vertex_descriptor V;

	vector<E> weak;
	Eit eit, elast;
	for (tie(eit, elast) = edges(g); eit != elast; ++eit) {
		E u1v2 = *eit;
		V u1 = source(u1v2, g), v2 = target(u1v2, g);
		if (out_degree(u1, g) != 2 || in_degree(v2, g) != 2)
			continue;

		Oit oit, olast;
		tie(oit, olast) = out_edges(u1, g);
		E u1v1;
		if (target(*oit, g) == v2) {
			++oit;
			u1v1 = *oit;
		} else {
			u1v1 = *oit;
			++oit;
		}
		assert(++oit == olast);
		V v1 = target(u1v1, g);
		assert(v1 != v2);
		if (in_degree(v1, g) != 1)
			continue;

		Iit iit, ilast;
		tie(iit, ilast) = in_edges(v2, g);
		E u2v2;
		if (source(*iit, g) == u1) {
			++iit;
			assert(iit != ilast);
			u2v2 = *iit;
		} else {
			assert(iit != ilast);
			u2v2 = *iit;
			++iit;
		}
		assert(++iit == ilast);
		V u2 = source(u2v2, g);
		assert(u1 != u2);
		if (out_degree(u2, g) != 1)
			continue;

		unsigned n = g[u1v2].numPairs;
		if (n < g[u1v1].numPairs && n < g[u2v2].numPairs)
			weak.push_back(u1v2);
	}

	if (opt::verbose > 1) {
		cerr << "Weak edges:\n";
		for (vector<E>::const_iterator it = weak.begin();
				it != weak.end(); ++it) {
			E e = *it;
			cerr << "\t\"" << source(e, g) << "\"->\"" << target(e, g)
				<< "\" [" << g[e] << "]\n";
		}
	}

	/** Remove the weak edges. */
	remove_edges(g, weak.begin(), weak.end());
	if (opt::verbose > 0) {
		cerr << "Removed " << weak.size() << " weak edges.\n";
		printGraphStats(cerr, g);
	}
}

/** Return whether the specified distance estimate is an exact
 * overlap.
 */
static bool isOverlap(const DistanceEst& d)
{
	if (d.stdDev == 0) {
		assert(d.distance < 0);
		return true;
	} else
		return false;
}

/** Add distance estimates to a path.
 * @param g0 the original graph
 * @param g1 the transformed graph
 */
static ContigPath addDistEst(const Graph& g0, const Graph& g1,
		const ContigPath& path)
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
		pair<E, bool> e0 = edge(u, v, g0);
		pair<E, bool> e1 = edge(u, v, g1);
		assert(e0.second || e1.second);
		const EP& ep = e0.second ? g0[e0.first] : g1[e1.first];
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

/** Read a graph from the specified file. */
static void readGraph(const string& path, Graph& g)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << path << "'...\n";
	ifstream fin(path.c_str());
	istream& in = path == "-" ? cin : fin;
	assert_good(in, path);
	read_graph(in, g, BetterDistanceEst());
	assert(in.eof());
	if (opt::verbose > 0)
		printGraphStats(cerr, g);
	ContigID::lock();
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
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	Graph g;
	if (optind < argc) {
		for (; optind < argc; optind++)
			readGraph(argv[optind], g);
	} else
		readGraph("-", g);

	// Add any missing complementary edges.
	addComplementaryEdges(g);
	Graph g0(g);

	// Filter the graph.
	filterGraph(g);
	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	// Resolve forks.
	resolveForks(g, g0);

	// Prune tips.
	pruneTips(g);

	// Remove repeats.
	removeRepeats(g);

	// Remove transitive edges.
	unsigned numTransitive = remove_transitive_edges(g);
	if (opt::verbose > 0) {
		cerr << "Removed " << numTransitive << " transitive edges.\n";
		printGraphStats(cerr, g);
	}

	// Prune tips.
	pruneTips(g);

	// Pop bubbles.
	vector<ContigNode> popped = popBubbles(g);
	if (opt::verbose > 0) {
		cerr << "Removed " << popped.size()
			<< " vertices in bubbles.\n";
		printGraphStats(cerr, g);
	}
	if (opt::verbose > 1) {
		cerr << "Popped: ";
		copy(popped.begin(), popped.end(),
				ostream_iterator<ContigNode>(cerr, " "));
		cerr << '\n';
	}

	// Remove weak edges.
	removeWeakEdges(g);
	Graph g1(g);

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
	ContigID::unlock();
	for (vector<ContigPath>::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		out << ContigID::create() << '\t'
			<< addDistEst(g0, g1, *it) << '\n';
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
