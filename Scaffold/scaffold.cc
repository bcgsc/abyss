#include "config.h"
#include "ContigNode.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "Estimate.h"
#include "IOUtil.h"
#include "Iterator.h"
#include "Uncompress.h"
#include "Graph/Assemble.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphAlgorithms.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include "Graph/PopBubbles.h"
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <utility>
#if _SQL
#include "DataBase/Options.h"
#include "DataBase/DB.h"
#endif

using namespace std;
using namespace std::rel_ops;
using boost::edge_bundle_type;
using boost::tie;

#define PROGRAM "abyss-scaffold"

#if _SQL
DB db;
#endif

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k<kmer> [OPTION]... FASTA|OVERLAP DIST...\n"
"Scaffold contigs using the distance estimate graph.\n"
"\n"
" Arguments:\n"
"\n"
"  FASTA    contigs in FASTA format\n"
"  OVERLAP  the contig overlap graph\n"
"  DIST     estimates of the distance between contigs\n"
"\n"
" Options:\n"
"\n"
"  -n, --npairs=N        minimum number of pairs [0]\n"
"  -s, --seed-length=N   minimum contig length [200]\n"
"          or -s N0-N1   Find the value of s in [N0,N1]\n"
"                        that maximizes the scaffold N50.\n"
"  -k, --kmer=N          length of a k-mer\n"
"      --min-gap=N       minimum scaffold gap length to output [50]\n"
"      --max-gap=N       maximum scaffold gap length to output [inf]\n"
"      --complex         remove complex transitive edges\n"
"      --no-complex      don't remove complex transitive edges [default]\n"
"      --SS              expect contigs to be oriented correctly\n"
"      --no-SS           no assumption about contig orientation [default]\n"
"  -o, --out=FILE        write the paths to FILE\n"
"  -g, --graph=FILE      write the graph to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
#if _SQL
"      --db=FILE         specify path of database repository in FILE\n"
"      --library=NAME    specify library NAME for sqlite\n"
"      --strain=NAME     specify strain NAME for sqlite\n"
"      --species=NAME    specify species NAME for sqlite\n"
#endif
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
#if _SQL
	string url;
	dbVars metaVars;
#endif
	unsigned k; // used by ContigProperties

	/** Minimum number of pairs. */
	static unsigned minNumPairs;

	/** Minimum contig length. */
	static unsigned minContigLength = 200;
	static unsigned minContigLengthEnd;

	/** Minimum scaffold gap length to output. */
	static int minGap = 50;

	/** Maximum scaffold gap length to output.
	 * -ve value means no maximum. */
	static int maxGap = -1;

	/** Write the paths to this file. */
	static string out;

	/** Write the graph to this file. */
	static string graphPath;

	/** Run a strand-specific RNA-Seq assembly. */
	static int ss;

	/** Verbose output. */
	int verbose; // used by PopBubbles

	/** Output format */
	int format = DOT; // used by DistanceEst

	/** Remove complex transitive edges */
	static int comp_trans;
}

static const char shortopts[] = "g:k:n:o:s:v";

#if _SQL
enum { OPT_HELP = 1, OPT_VERSION, OPT_MIN_GAP, OPT_MAX_GAP, OPT_COMP,
	OPT_DB, OPT_LIBRARY, OPT_STRAIN, OPT_SPECIES };
#else
enum { OPT_HELP = 1, OPT_VERSION, OPT_MIN_GAP, OPT_MAX_GAP, OPT_COMP };
#endif

static const struct option longopts[] = {
	{ "graph",       no_argument,       NULL, 'g' },
	{ "kmer",        required_argument, NULL, 'k' },
	{ "min-gap",     required_argument, NULL, OPT_MIN_GAP },
	{ "max-gap",     required_argument, NULL, OPT_MAX_GAP },
	{ "npairs",      required_argument, NULL, 'n' },
	{ "out",         required_argument, NULL, 'o' },
	{ "seed-length", required_argument, NULL, 's' },
	{ "complex",     no_argument, &opt::comp_trans, 1 },
	{ "no-complex",  no_argument, &opt::comp_trans, 0 },
	{ "SS",          no_argument,       &opt::ss, 1 },
	{ "no-SS",       no_argument,       &opt::ss, 0 },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
#if _SQL
	{ "db",          required_argument, NULL, OPT_DB },
	{ "library",     required_argument, NULL, OPT_LIBRARY },
	{ "strain",      required_argument, NULL, OPT_STRAIN },
	{ "species",     required_argument, NULL, OPT_SPECIES },
#endif
	{ NULL, 0, NULL, 0 }
};

/** A distance estimate graph. */
typedef DirectedGraph<Length, DistanceEst> DG;
typedef ContigGraph<DG> Graph;

/** Return whether this edge is invalid.
 * An edge is invalid when the overlap is larger than the length of
 * either of its incident sequences.
 */
struct InvalidEdge {
	InvalidEdge(Graph& g) : m_g(g) { }
	bool operator()(graph_traits<Graph>::edge_descriptor e) const
	{
		int d = m_g[e].distance;
		int ulength = m_g[source(e, m_g)].length;
		int vlength = m_g[target(e, m_g)].length;
		return d + ulength <= 0 || d + vlength <= 0;
	}
	const Graph& m_g;
};

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
static void filterGraph(Graph& g, unsigned minContigLength)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_descriptor V;
	typedef GTraits::vertex_iterator Vit;

	// Remove short contigs.
	unsigned numRemovedV = 0;
	std::pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		if (g[u].length < minContigLength)
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
#if _SQL
	addToDb(db, "V_removed", numRemovedV);
	addToDb(db, "E_removed", numRemovedE);
#endif
}

/** Return true if the specified edge is a cycle. */
static bool isCycle(Graph& g, graph_traits<Graph>::edge_descriptor e)
{
	return edge(target(e, g), source(e, g), g).second;
}

/** Remove simple cycles of length two from the graph. */
static void removeCycles(Graph& g)
{
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef graph_traits<Graph>::edge_iterator Eit;

	// Identify the cycles.
	vector<E> cycles;
	Eit eit, elast;
	for (tie(eit, elast) = edges(g); eit != elast; ++eit) {
		E e = *eit;
		if (isCycle(g, e))
			cycles.push_back(e);
	}

	/** Remove the cycles. */
	remove_edges(g, cycles.begin(), cycles.end());
	if (opt::verbose > 0) {
		cerr << "Removed " << cycles.size() << " cyclic edges.\n";
		printGraphStats(cerr, g);
	}
#if _SQL
	addToDb(db, "E_removed_cyclic", cycles.size());
#endif
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
						cerr << "cycle: " << get(vertex_name, g, v1)
							<< ' ' << get(vertex_name, g, v2) << '\n';
				} else if (e12.second || e21.second) {
					E e = e12.second ? e12.first : e21.first;
					V v = source(e, g0), w = target(e, g0);
					add_edge(v, w, g0[e], g);
					numEdges++;
					if (opt::verbose > 1)
						cerr << get(vertex_name, g, u)
							<< " -> " << get(vertex_name, g, v)
							<< " -> " << get(vertex_name, g, w)
							<< " [" << g0[e] << "]\n";
				}
			}
		}
	}
	if (opt::verbose > 0)
		cerr << "Added " << numEdges
			<< " edges to ambiguous vertices.\n";
#if _SQL
	addToDb(db, "E_added_ambig", numEdges);
#endif
}

/** Remove tips.
 * For an edge (u,v), remove the vertex v if deg+(u) > 1
 * and deg-(v) = 1 and deg+(v) = 0.
 */
static void pruneTips(Graph& g)
{
	/** Identify the tips. */
	size_t n = 0;
	pruneTips(g, CountingOutputIterator(n));

	if (opt::verbose > 0) {
		cerr << "Removed " << n << " tips.\n";
		printGraphStats(cerr, g);
	}
#if _SQL
	addToDb(db, "Tips_removed", n);
#endif
}

/** Remove repetitive vertices from this graph.
 * input: digraph g { t1->v1 t2->v2 t1->u t2->u u->v1 u->v2 }
 * operation: remove vertex u
 * output: digraph g { t1->v1 t2->v2 }
 */
static void removeRepeats(Graph& g)
{
	typedef graph_traits<Graph>::adjacency_iterator Ait;
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef graph_traits<Graph>::vertex_descriptor V;

	vector<V> repeats;
	vector<E> transitive;
	find_transitive_edges(g, back_inserter(transitive));
	for (vector<E>::const_iterator it = transitive.begin();
			it != transitive.end(); ++it) {
		// Iterate through the transitive edges, u->w1.
		V u = source(*it, g), w1 = target(*it, g);
		Ait vit, vlast;
		for (tie(vit, vlast) = adjacent_vertices(u, g);
				vit != vlast; ++vit) {
			V v = *vit;
			assert(u != v); // no self loops
			if (!edge(v, w1, g).second)
				continue;
			// u->w1 is a transitive edge spanning u->v->w1.
			Ait wit, wlast;
			for (tie(wit, wlast) = adjacent_vertices(v, g);
					wit != wlast; ++wit) {
				// For each edge v->w2, check that an edge
				// w1->w2 or w2->w1 exists. If not, v is a repeat.
				V w2 = *wit;
				assert(v != w2); // no self loops
				if (w1 != w2
						&& !edge(w1, w2, g).second
						&& !edge(w2, w1, g).second) {
					repeats.push_back(v);
					break;
				}
			}
		}
	}

	sort(repeats.begin(), repeats.end());
	repeats.erase(unique(repeats.begin(), repeats.end()),
			repeats.end());
	if (opt::verbose > 1) {
		cerr << "Ambiguous:";
		for (vector<V>::const_iterator it = repeats.begin();
				it != repeats.end(); ++it)
			cerr << ' ' << get(vertex_name, g, *it);
		cerr << '\n';
	}

	// Remove the repetitive vertices.
	unsigned numRemoved = 0;
	for (vector<V>::const_iterator it = repeats.begin();
			it != repeats.end(); ++it) {
		V u = *it;
		V uc = get(vertex_complement, g, u);
		clear_out_edges(u, g);
		if (it != repeats.begin() && it[-1] == uc) {
			remove_vertex(u, g);
			numRemoved++;
		}
	}

	if (opt::verbose > 0) {
		cerr << "Cleared "
			<< repeats.size() << " ambiguous vertices.\n"
			<< "Removed "
			<< numRemoved << " ambiguous vertices.\n";
		printGraphStats(cerr, g);
	}
#if _SQL
	addToDb(db, "V_cleared_ambg", repeats.size());
	addToDb(db, "V_removed_ambg", numRemoved);
#endif
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
			cerr << '\t' << get(edge_name, g, e)
				<< " [" << g[e] << "]\n";
		}
	}

	/** Remove the weak edges. */
	remove_edges(g, weak.begin(), weak.end());
	if (opt::verbose > 0) {
		cerr << "Removed " << weak.size() << " weak edges.\n";
		printGraphStats(cerr, g);
	}
#if _SQL
	addToDb(db, "E_removed_weak", weak.size());
#endif
}

static void removeLongEdges(Graph& g)
{
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef graph_traits<Graph>::edge_iterator Eit;

	vector<E> long_e;
	Eit eit, elast;
	for (tie(eit, elast) = edges(g); eit != elast; ++eit) {
		E e = *eit;
		if (g[e].distance > opt::maxGap)
			long_e.push_back(e);
	}
	remove_edges(g, long_e.begin(), long_e.end());
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

#if _SQL
	vector<int> vals = passGraphStatsVal(g);
	vector<string> keys = make_vector<string>()
		<< "V_readGraph"
		<< "E_readGraph"
		<< "degree0_readGraph"
		<< "degree1_readGraph"
		<< "degree234_readGraph"
		<< "degree5_readGraph"
		<< "max_readGraph";

	for(unsigned i=0; i<vals.size(); i++)
		addToDb(db, keys[i], vals[i]);
#endif

	g_contigNames.lock();
}

/** Return the scaffold length of [first, last), not counting gaps. */
template<typename It>
unsigned addLength(const Graph& g, It first, It last)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	assert(first != last);
	unsigned length = g[*first].length;
	for (It it = first + 1; it != last; ++it) {
		V u = *(it - 1);
		V v = *it;
		length += min(0, get(edge_bundle, g, u, v).distance);
		length += g[v].length;
	}
	return length;
}

/** A container of contig paths. */
typedef vector<ContigPath> ContigPaths;

/** Build the scaffold length histogram. */
static Histogram buildScaffoldLengthHistogram(
		Graph& g, const ContigPaths& paths)
{
	Histogram h;

	// Clear the removed flag.
	typedef graph_traits<Graph>::vertex_iterator Vit;
	Vit uit, ulast;
	for (tie(uit, ulast) = vertices(g); uit != ulast; ++uit)
		put(vertex_removed, g, *uit, false);

	// Remove the vertices that are used in paths
	// and add the lengths of the scaffolds.
	for (ContigPaths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		h.insert(addLength(g, it->begin(), it->end()));
		remove_vertex_if(g, it->begin(), it->end(),
				not1(std::mem_fun_ref(&ContigNode::ambiguous)));
	}

	// Add the contigs that were not used in paths.
	for (tie(uit, ulast) = vertices(g); uit != ulast; ++++uit) {
		typedef graph_traits<Graph>::vertex_descriptor V;
		V u = *uit;
		if (!get(vertex_removed, g, u))
			h.insert(g[u].length);
	}

	return h;
}

#if _SQL
/** Add contiguity stats to database */
static void addCntgStatsToDb(
		const Histogram h, const unsigned min)
{
	vector<int> vals = passContiguityStatsVal(h, min);
	vector<string> keys = make_vector<string>()
		<< "n"
		<< "n200"
		<< "nN50"
		<< "min"
		<< "N80"
		<< "N50"
		<< "N20"
		<< "Esize"
		<< "max"
		<< "sum"
		<< "nNG50"
		<< "NG50";

	for(unsigned i=0; i<vals.size(); i++)
		addToDb(db, keys[i], vals[i]);
}
#endif

/** Build scaffold paths.
 * @param output write the results
 * @return the scaffold N50
 */
unsigned scaffold(const Graph& g0, unsigned minContigLength,
		bool output)
{
	Graph g(g0);

	// Filter the graph.
	filterGraph(g, minContigLength);
	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	// Remove cycles.
	removeCycles(g);

	// Resolve forks.
	resolveForks(g, g0);

	// Prune tips.
	pruneTips(g);

	// Remove repeats.
	removeRepeats(g);

	// Remove transitive edges.
	unsigned numTransitive;
	if (opt::comp_trans)
		numTransitive = remove_complex_transitive_edges(g);
	else
		numTransitive = remove_transitive_edges(g);

	if (opt::verbose > 0) {
		cerr << "Removed " << numTransitive << " transitive edges.\n";
		printGraphStats(cerr, g);
	}

#if _SQL
	addToDb(db, "Edges_transitive", numTransitive);
#endif
	// Prune tips.
	pruneTips(g);

	// Pop bubbles.
	typedef graph_traits<Graph>::vertex_descriptor V;
	vector<V> popped = popBubbles(g);
	if (opt::verbose > 0) {
		cerr << "Removed " << popped.size()
			<< " vertices in bubbles.\n";
		printGraphStats(cerr, g);
	}
#if _SQL
	addToDb(db, "Vertices_bubblePopped", popped.size());
#endif
	if (opt::verbose > 1) {
		cerr << "Popped:";
		for (vector<V>::const_iterator it = popped.begin();
				it != popped.end(); ++it)
			cerr << ' ' << get(vertex_name, g, *it);
		cerr << '\n';
	}

	// Remove weak edges.
	removeWeakEdges(g);

	// Remove any edges longer than opt::maxGap.
	if (opt::maxGap >= 0)
		removeLongEdges(g);

	// Assemble the paths.
	ContigPaths paths;
	assembleDFS(g, back_inserter(paths), opt::ss);
	sort(paths.begin(), paths.end());
	unsigned n = 0;
	if (opt::verbose > 0) {
		for (ContigPaths::const_iterator it = paths.begin();
				it != paths.end(); ++it)
			n += it->size();
		cerr << "Assembled " << n << " contigs in "
			<< paths.size() << " scaffolds.\n";
		printGraphStats(cerr, g);
	}
#if _SQL
	addToDb(db, "contigs_assembled", n);
	addToDb(db, "scaffolds_assembled", paths.size());
#endif

	const unsigned STATS_MIN_LENGTH = opt::minContigLength;
	if (!output) {
		static bool printHeader = true;
		Histogram h = buildScaffoldLengthHistogram(g, paths);
		printContiguityStats(cerr, h, STATS_MIN_LENGTH,
				printHeader)
			<< "\ts=" << minContigLength << '\n';
		if (opt::verbose == 0)
			printHeader = false;
#if _SQL
		addCntgStatsToDb(h, STATS_MIN_LENGTH);
#endif
		return h.trimLow(STATS_MIN_LENGTH).n50();
	}

	// Output the paths.
	ofstream fout(opt::out.c_str());
	ostream& out = opt::out.empty() || opt::out == "-" ? cout : fout;
	assert_good(out, opt::out);
	g_contigNames.unlock();
	for (vector<ContigPath>::const_iterator it = paths.begin();
			it != paths.end(); ++it)
		out << createContigName() << '\t'
			<< addDistEst(g0, g, *it) << '\n';
	assert_good(out, opt::out);

	// Output the graph.
	if (!opt::graphPath.empty()) {
		ofstream out(opt::graphPath.c_str());
		assert_good(out, opt::graphPath);
		write_dot(out, g);
		assert_good(out, opt::graphPath);
	}

	// Print assembly contiguity statistics.
	Histogram h = buildScaffoldLengthHistogram(g, paths);
	printContiguityStats(cerr, h, STATS_MIN_LENGTH) << '\n';
#if _SQL
	addCntgStatsToDb(h, STATS_MIN_LENGTH);
#endif
	return h.trimLow(STATS_MIN_LENGTH).n50();
}

/** Run abyss-scaffold. */
int main(int argc, char** argv)
{
#if _SQL
	opt::metaVars.resize(3);
#endif

	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			die = true;
			break;
		  case 'k':
			arg >> opt::k;
			break;
		  case 'g':
			arg >> opt::graphPath;
			break;
		  case 'n':
			arg >> opt::minNumPairs;
			break;
		  case 'o':
			arg >> opt::out;
			break;
		  case 's':
			arg >> opt::minContigLength;
			if (arg.peek() == '-') {
				opt::minContigLengthEnd = 100 * opt::minContigLength;
				arg >> expect("-") >> opt::minContigLengthEnd;
				assert(opt::minContigLength
						<= opt::minContigLengthEnd);
			}
			break;
		  case 'v':
			opt::verbose++;
			break;
		  case OPT_MIN_GAP:
			arg >> opt::minGap;
			break;
		  case OPT_MAX_GAP:
			arg >> opt::maxGap;
			break;
		  case OPT_HELP:
			cout << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case OPT_VERSION:
			cout << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
#if _SQL
		  case OPT_DB:
			arg >> opt::url;
			break;
		  case OPT_LIBRARY:
			arg >> opt::metaVars[0];
			break;
		  case OPT_STRAIN:
			arg >> opt::metaVars[1];
			break;
		  case OPT_SPECIES:
			arg >> opt::metaVars[2];
			break;
#endif
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
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
#if _SQL
	init(db,
			opt::url,
			opt::verbose,
			PROGRAM,
			opt::getCommand(argc, argv),
			opt::metaVars
	);
	addToDb(db, "K", opt::k);
#endif

	Graph g;
	if (optind < argc) {
		for (; optind < argc; optind++)
			readGraph(argv[optind], g);
	} else
		readGraph("-", g);

	// Add any missing complementary edges.
	size_t numAdded = addComplementaryEdges(g);
	if (opt::verbose > 0) {
		cerr << "Added " << numAdded << " complementary edges.\n";
		printGraphStats(cerr, g);
	}

#if _SQL
	addToDb(db, "add_complement_edges", numAdded);
#endif

	// Remove invalid edges.
	unsigned numBefore = num_edges(g);
	remove_edge_if(InvalidEdge(g), static_cast<DG&>(g));
	unsigned numRemoved = numBefore - num_edges(g);
	if (numRemoved > 0)
		cerr << "warning: Removed "
			<< numRemoved << " invalid edges.\n";

#if _SQL
	addToDb(db, "Edges_invalid", numRemoved);
#endif
	if (opt::minContigLengthEnd == 0) {
		scaffold(g, opt::minContigLength, true);
		return 0;
	}

	// Find the value of s that maximizes the scaffold N50.
	unsigned bests = 0, bestN50 = 0;
	const double STEP = cbrt(10); // Three steps per decade.
	unsigned ilast = (unsigned)round(
			log(opt::minContigLengthEnd) / log(STEP));
	for (unsigned i = (unsigned)round(
				log(opt::minContigLength) / log(STEP));
			i <= ilast; ++i) {
		unsigned s = (unsigned)pow(STEP, (int)i);

		// Round to 1 figure.
		double nearestDecade = pow(10, floor(log10(s)));
		s = unsigned(round(s / nearestDecade) * nearestDecade);

		unsigned n50 = scaffold(g, s, false);
		if (opt::verbose > 0)
			cerr << '\n';
		if (n50 > bestN50) {
			bestN50 = n50;
			bests = s;
		}
	}

	bestN50 = scaffold(g, bests, true);
	cerr << "Best scaffold N50 is " << bestN50
		<< " at s=" << bests << ".\n";

	return 0;
}
