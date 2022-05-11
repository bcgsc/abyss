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
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include <algorithm>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
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
using namespace boost::lambda;
using boost::tie;

#define PROGRAM "abyss-filtergraph"

static const char VERSION_MESSAGE[] =
    PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
            "Written by Tony Raymond.\n"
            "\n"
            "Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " -k<kmer> [OPTION]... ADJ [FASTA]\n"
    "Remove short contigs that do not contribute any relevant\n"
    "information to the assembly.\n"
    "\n"
    " Arguments:\n"
    "\n"
    "  ADJ    contig adjacency graph\n"
    "  FASTA  contigs to check consistency of ADJ edges\n"
    "\n"
    " Options:\n"
    "\n"
    "  -k, --kmer=N            k-mer size\n"
    "      --SS                expect contigs to be oriented correctly\n"
    "      --no-SS             no assumption about contig orientation\n"
    "  -T, --island=N          remove islands shorter than N [0]\n"
    "  -t, --tip=N             remove tips shorter than N [0]\n"
    "  -l, --length=N          remove contigs shorter than N [0]\n"
    "  -L, --max-length=N      remove contigs longer than N [0]\n"
    "  -c, --coverage=FLOAT    remove contigs with mean k-mer coverage less than FLOAT [0]\n"
    "  -C, --max-coverage=FLOAT remove contigs with mean k-mer coverage at least FLOAT [0]\n"
    "      --shim              remove filler contigs that only contribute\n"
    "                          to adjacency [default]\n"
    "      --no-shim           disable filler contigs removal\n"
    "      --shim-max-degree=N only remove shims where the smaller of \n"
    "                          in/out degree is smaller than N [1]\n"
    "  -m, --min-overlap=N     require a minimum overlap of N bases [10]\n"
    "      --assemble          assemble unambiguous paths\n"
    "      --no-assemble       disable assembling of paths [default]\n"
    "  -g, --graph=FILE        write the contig adjacency graph to FILE\n"
    "  -i, --ignore=FILE       ignore contigs seen in FILE\n"
    "  -r, --remove=FILE       remove contigs seen in FILE\n"
    "      --adj               output the graph in ADJ format [default]\n"
    "      --asqg              output the graph in ASQG format\n"
    "      --dot               output the graph in GraphViz format\n"
    "      --gfa               output the graph in GFA1 format\n"
    "      --gfa1              output the graph in GFA1 format\n"
    "      --gfa2              output the graph in GFA2 format\n"
    "      --gv                output the graph in GraphViz format\n"
    "      --sam               output the graph in SAM format\n"
    "  -v, --verbose           display verbose output\n"
    "      --help              display this help and exit\n"
    "      --version           output version information and exit\n"
    "\n"
    "Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
unsigned k; // used by ContigProperties

/** Run a strand-specific RNA-Seq assembly. */
static int ss;

/** Remove island contigs less than this length. */
static unsigned minIslandLen = 0;

/** Remove tips less than this length. */
static unsigned minTipLen = 0;

/** Remove all contigs less than this length. */
static unsigned minLen = 0;

/** Remove all contigs more than this length. */
static unsigned maxLen = 0;

/** Remove contigs with mean k-mer coverage less than this threshold. */
static float minCoverage = 0;

/** Remove contigs with mean k-mer coverage at least this threshold. */
static float maxCoverage = 0;

/** Remove short contigs that don't contribute any sequence. */
static int shim = 1;

/** Only remove shims where the smaller of in/out degree is small
 * enough. */
static unsigned shimMaxDegree = 1;

/** Assemble unambiguous paths. */
static int assemble = 0;

/** Write the contig adjacency graph to this file. */
static string graphPath;

/** Contigs to ignore. */
static string ignorePath;

/** Contigs to remove. */
static string removePath;

/** The minimum overlap allowed between two contigs. */
static int minOverlap = 10;

/** Output graph format. */
int format = ADJ; // used by ContigProperties
}

static const char shortopts[] = "c:C:g:i:r:k:l:L:m:t:T:v";

enum
{
	OPT_HELP = 1,
	OPT_VERSION,
	OPT_SHIM_MAX_DEG
};

static const struct option longopts[] = {
	{ "adj", no_argument, &opt::format, ADJ },
	{ "asqg", no_argument, &opt::format, ASQG },
	{ "dot", no_argument, &opt::format, DOT },
	{ "gfa", no_argument, &opt::format, GFA1 },
	{ "gfa1", no_argument, &opt::format, GFA1 },
	{ "gfa2", no_argument, &opt::format, GFA2 },
	{ "gv", no_argument, &opt::format, DOT },
	{ "sam", no_argument, &opt::format, SAM },
	{ "graph", required_argument, NULL, 'g' },
	{ "ignore", required_argument, NULL, 'i' },
	{ "remove", required_argument, NULL, 'r' },
	{ "SS", no_argument, &opt::ss, 1 },
	{ "no-SS", no_argument, &opt::ss, 0 },
	{ "kmer", required_argument, NULL, 'k' },
	{ "island", required_argument, NULL, 'T' },
	{ "tip", required_argument, NULL, 't' },
	{ "length", required_argument, NULL, 'l' },
	{ "max-length", required_argument, NULL, 'L' },
	{ "coverage", required_argument, NULL, 'c' },
	{ "max-coverage", required_argument, NULL, 'C' },
	{ "shim", no_argument, &opt::shim, 1 },
	{ "no-shim", no_argument, &opt::shim, 0 },
	{ "shim-max-degree", required_argument, NULL, OPT_SHIM_MAX_DEG },
	{ "assemble", no_argument, &opt::assemble, 1 },
	{ "no-assemble", no_argument, &opt::assemble, 0 },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "verbose", no_argument, NULL, 'v' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static vector<ContigID> g_removed;

/** Contig adjacency graph. */
typedef ContigGraph<DirectedGraph<ContigProperties, Distance>> Graph;
typedef Graph::vertex_descriptor vertex_descriptor;
typedef Graph::edge_descriptor edge_descriptor;

/** Data for verbose output. */
static struct
{
	unsigned removed;
	unsigned tails;
	unsigned too_long;
	unsigned too_complex;
	unsigned self_adj;
	unsigned checked;
	unsigned parallel_edge;
} g_count;

/** Returns if the contig can be removed from the graph. */
static bool
removable(const Graph* pg, vertex_descriptor v)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::out_edge_iterator OEit;
	typedef GTraits::in_edge_iterator IEit;
	typedef GTraits::vertex_descriptor V;

	const Graph& g = *pg;

	// Check if previously removed
	if (get(vertex_removed, g, v)) {
		g_count.removed++;
		return false;
	}

	unsigned min_degree = min(out_degree(v, g), in_degree(v, g));

	// Check for tails
	if (min_degree == 0) {
		g_count.tails++;
		return false;
	}

	// Check that the result will be less complex that the original
	if (min_degree > opt::shimMaxDegree) {
		g_count.too_complex++;
		return false;
	}

	// Check if self adjacent
	OEit oei0, oei1;
	tie(oei0, oei1) = out_edges(v, g);
	for (OEit vw = oei0; vw != oei1; ++vw) {
		V w = target(*vw, g);
		V vc = get(vertex_complement, g, v);
		if (v == w || vc == w) {
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

	if (g[*maxuv].distance + (int)g[v].length + g[*maxvw].distance > -opt::minOverlap) {
		g_count.too_long++;
		return false;
	}
	return true;
}

/** Data to store information of an edge. */
struct EdgeInfo
{
	vertex_descriptor u;
	vertex_descriptor w;
	edge_bundle_type<Graph>::type ep;

	EdgeInfo(vertex_descriptor u, vertex_descriptor w, int ep)
	  : u(u)
	  , w(w)
	  , ep(ep)
	{}
	EdgeInfo()
	  : u()
	  , w()
	  , ep()
	{}
};

/** Returns a list of edges that may be added when the vertex v is
 * removed. */
static bool
findNewEdges(
    const Graph& g,
    vertex_descriptor v,
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
			int x = g[*uv].distance + (int)g[v].length + g[*vw].distance;
			assert(x <= 0);
			EdgeInfo ed(source(*uv, g), target(*vw, g), x);
			eds.push_back(ed);
			if (out_degree(v, g) > 1)
				marked.push_back(ed.u);
			if (in_degree(v, g) > 1)
				marked.push_back(ed.w);
		}
	}
	for (vector<V>::const_iterator it = marked.begin(); it != marked.end(); it++)
		markedContigs[get(vertex_index, g, *it)] = true;
	return true;
}

/** Adds all edges described in the vector eds. */
static void
addNewEdges(Graph& g, const vector<EdgeInfo>& eds)
{
	for (vector<EdgeInfo>::const_iterator edsit = eds.begin(); edsit != eds.end(); ++edsit) {
		// Don't add parallel edges! This can happen when removing a palindrome.
		if (edge(edsit->u, edsit->w, g).second) {
			g_count.parallel_edge++;
			continue;
		}
		assert(edsit->ep.distance <= -opt::minOverlap);
		add_edge(edsit->u, edsit->w, edsit->ep, g);
	}
}

static void
removeContig(vertex_descriptor v, Graph& g)
{
	clear_vertex(v, g);
	remove_vertex(v, g);
	g_removed.push_back(get(vertex_contig_index, g, v));
	g_count.removed++;
}

/** Remove the specified contig from the adjacency graph. */
static void
removeContigs(Graph& g, vector<vertex_descriptor>& sc)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_descriptor V;

	vector<vertex_descriptor> out;
	out.reserve(sc.size());

	vector<bool> markedContigs(g.num_vertices());
	for (vector<vertex_descriptor>::iterator it = sc.begin(); it != sc.end(); ++it) {
		V v = *it;
		if (opt::verbose > 0 && ++g_count.checked % 10000000 == 0)
			cerr << "Removed " << g_count.removed << "/" << g_count.checked
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

		removeContig(v, g);
	}
	sc.swap(out);
}

/** Return the value of the bit at the specified index. */
struct Marked
{
	typedef vector<bool> Data;
	Marked(const Graph& g, const Data& data)
	  : m_g(g)
	  , m_data(data)
	{}
	bool operator()(vertex_descriptor u) const { return m_data[get(vertex_contig_index, m_g, u)]; }

	typedef vertex_descriptor argument_type;
	typedef bool result_type;

  private:
	const Graph& m_g;
	const Data& m_data;
};

/** Finds all potentially removable contigs in the graph. */
static void
findShortContigs(const Graph& g, const vector<bool>& seen, vector<vertex_descriptor>& sc)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_iterator Vit;
	Vit first, second;
	tie(first, second) = vertices(g);
	::copy_if(
	    first,
	    second,
	    back_inserter(sc),
	    !boost::lambda::bind(Marked(g, seen), _1) && boost::lambda::bind(removable, &g, _1));
}

/** Functor used for sorting contigs based on degree, then size,
 * and then ID. */
struct sortContigs
{
	const Graph& g;

	sortContigs(const Graph& g)
	  : g(g)
	{}

	template<typename V>
	bool operator()(V a, V b)
	{
		const ContigProperties& ap = g[a];
		const ContigProperties& bp = g[b];

		unsigned dega = out_degree(a, g) * in_degree(a, g);
		unsigned degb = out_degree(b, g) * in_degree(b, g);

		return dega != degb ? dega < degb : ap.length != bp.length ? ap.length < bp.length : a < b;
	}
};

struct ShorterThanX
{
	const Graph& g;
	const vector<bool>& seen;
	size_t x;

	ShorterThanX(const Graph& g, const vector<bool>& seen, size_t x)
	  : g(g)
	  , seen(seen)
	  , x(x)
	{}

	bool operator()(vertex_descriptor y) const
	{
		return g[y].length < x && !get(vertex_removed, g, y) &&
		       !seen[get(vertex_contig_index, g, y)];
	}

	typedef vertex_descriptor argument_type;
	typedef bool result_type;
};

struct LongerThanX
{
	const Graph& g;
	const vector<bool>& seen;
	size_t x;

	LongerThanX(const Graph& g, const vector<bool>& seen, size_t x)
	  : g(g)
	  , seen(seen)
	  , x(x)
	{}

	bool operator()(vertex_descriptor y) const
	{
		return g[y].length > x && !get(vertex_removed, g, y) &&
		       !seen[get(vertex_contig_index, g, y)];
	}

	typedef vertex_descriptor argument_type;
	typedef bool result_type;
};

struct CoverageLessThan
{
	const Graph& g;
	const vector<bool>& seen;
	float minCov;

	CoverageLessThan(const Graph& g, const vector<bool>& seen, float minCov)
	  : g(g)
	  , seen(seen)
	  , minCov(minCov)
	{}

	bool operator()(vertex_descriptor u) const
	{
		assert(opt::k > 0);
		float meanCoverage = (float)g[u].coverage / (g[u].length - opt::k + 1);
		return meanCoverage < minCov && !get(vertex_removed, g, u) &&
		       !seen[get(vertex_contig_index, g, u)];
	}

	typedef vertex_descriptor argument_type;
	typedef bool result_type;
};

static void
removeShims(Graph& g, const vector<bool>& seen)
{
	if (opt::verbose > 0)
		cerr << "Removing shim contigs from the graph...\n";
	vector<vertex_descriptor> shortContigs;
	findShortContigs(g, seen, shortContigs);
	for (unsigned i = 0; !shortContigs.empty(); ++i) {
		if (opt::verbose > 0)
			cerr << "Pass " << i + 1 << ": Checking " << shortContigs.size() << " contigs.\n";
		sort(shortContigs.begin(), shortContigs.end(), sortContigs(g));
		removeContigs(g, shortContigs);
	}
	if (opt::verbose > 0) {
		cerr << "Shim removal stats:\n";
		cerr << "Removed: " << g_count.removed / 2 << " Too Complex: " << g_count.too_complex / 2
		     << " Tails: " << g_count.tails / 2 << " Too Long: " << g_count.too_long / 2
		     << " Self Adjacent: " << g_count.self_adj / 2
		     << " Parallel Edges: " << g_count.parallel_edge / 2 << '\n';
	}
}

template<typename pred>
static void
removeContigs_if(Graph& g, pred p)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::vertex_iterator Vit;
	typedef GTraits::vertex_descriptor V;
	Vit first, second;
	tie(first, second) = vertices(g);
	vector<V> sc;
	::copy_if(first, second, back_inserter(sc), p);
	remove_vertex_if(g, sc.begin(), sc.end(), True<V>());
	transform(sc.begin(), sc.end(), back_inserter(g_removed), [](const ContigNode& c) {
		return c.contigIndex();
	});
	if (opt::verbose > 0)
		cerr << "Removed " << sc.size() / 2 << " contigs.\n";
}

/** Contig sequences. */
typedef vector<const_string> Contigs;
static Contigs g_contigs;

/** Return the sequence of vertex u. */
static string
getSequence(const Graph& g, vertex_descriptor u)
{
	size_t i = get(vertex_contig_index, g, u);
	assert(i < g_contigs.size());
	string seq(g_contigs[i]);
	return get(vertex_sense, g, u) ? reverseComplement(seq) : seq;
}

/** Return whether the specified edge is inconsistent. */
struct is_edge_inconsistent
{
	const Graph& g;

	is_edge_inconsistent(const Graph& g)
	  : g(g)
	{}

	bool operator()(edge_descriptor e) const
	{
		vertex_descriptor u = source(e, g);
		vertex_descriptor v = target(e, g);

		int overlap = g[e].distance;
		assert(overlap < 0);

		string su = getSequence(g, u);
		string sv = getSequence(g, v);
		const unsigned u_start = su.length() + overlap;

		for (unsigned i = 0; i < (unsigned)-overlap; i++)
			if (!(ambiguityToBitmask(su[u_start + i]) & ambiguityToBitmask(sv[i])))
				return true;
		return false;
	}

	typedef edge_descriptor argument_type;
	typedef bool result_type;
};

template<typename It>
static void
remove_edge(Graph& g, It first, It last)
{
	for (; first != last; first++)
		remove_edge(*first, g);
}

template<typename pred>
static void
removeEdges_if(Graph& g, pred p)
{
	typedef graph_traits<Graph> GTraits;
	typedef GTraits::edge_iterator Eit;
	typedef GTraits::edge_descriptor E;
	Eit first, second;
	tie(first, second) = edges(g);
	vector<E> sc;
	::copy_if(first, second, back_inserter(sc), p);
	remove_edge(g, sc.begin(), sc.end());
	if (opt::verbose > 0) {
		cerr << "Edge removal stats:\n";
		cerr << "Removed: " << sc.size() << '\n';
	}
}

int
main(int argc, char** argv)
{
	string commandLine;
	{
		ostringstream ss;
		char** last = argv + argc - 1;
		copy(argv, last, ostream_iterator<const char*>(ss, " "));
		ss << *last;
		commandLine = ss.str();
	}

	bool die = false;
	for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case '?':
			die = true;
			break;
		case 'c':
			arg >> opt::minCoverage;
			break;
		case 'C':
			arg >> opt::maxCoverage;
			break;
		case 'l':
			arg >> opt::minLen;
			break;
		case 'L':
			arg >> opt::maxLen;
			break;
		case 'm':
			arg >> opt::minOverlap;
			break;
		case 'g':
			arg >> opt::graphPath;
			break;
		case 'i':
			arg >> opt::ignorePath;
			break;
		case 'r':
			arg >> opt::removePath;
			break;
		case 'k':
			arg >> opt::k;
			break;
		case 'T':
			arg >> opt::minIslandLen;
			break;
		case 't':
			arg >> opt::minTipLen;
			break;
		case 'v':
			opt::verbose++;
			break;
		case OPT_SHIM_MAX_DEG:
			arg >> opt::shimMaxDegree;
			break;
		case OPT_HELP:
			cout << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		case OPT_VERSION:
			cout << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-" << (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (opt::minOverlap < 0) {
		cerr << PROGRAM ": "
		     << "--min-overlap must be a positive integer.\n";
		die = true;
	}

	if (opt::k <= 0) {
		cerr << PROGRAM ": "
		     << "missing -k,--kmer option\n";
		die = true;
	}

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (argc - optind > 2) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	Graph g;
	// Read the contig adjacency graph.
	{
		string adjPath(argv[optind++]);
		if (opt::verbose > 0)
			cerr << "Loading graph from file: " << adjPath << '\n';
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

	// Remove list of contigs
	if (!opt::removePath.empty()) {
		ifstream in(opt::removePath.c_str());
		assert(in.good());
		string s;
		size_t b = g_removed.size();
		while (in >> s) {
			size_t i = get(g_contigNames, s);
			removeContig(ContigNode(i, 0), g);
		}
		assert(in.eof());
		if (opt::verbose)
			cerr << "Removed " << g_removed.size() - b << " contigs.\n";
	}

	// Remove shims.
	if (opt::shim)
		removeShims(g, seen);

	// Remove islands.
	if (opt::minIslandLen > 0) {
		size_t s = g_removed.size();
		removeIslands_if(g, back_inserter(g_removed), ShorterThanX(g, seen, opt::minIslandLen));
		if (opt::verbose)
			cerr << "Removed " << g_removed.size() - s << " islands.\n";
	}

	// Remove tips.
	if (opt::minTipLen > 0) {
		size_t s, prev;
		s = g_removed.size();
		do {
			prev = g_removed.size();
			pruneTips_if(g, back_inserter(g_removed), ShorterThanX(g, seen, opt::minTipLen));
		} while (prev < g_removed.size());
		if (opt::verbose)
			cerr << "Removed " << g_removed.size() - s << " tips.\n";
	}

	// Remove short contigs.
	if (opt::minLen > 0)
		removeContigs_if(g, ShorterThanX(g, seen, opt::minLen));

	// Remove long contigs.
	if (opt::maxLen > 0)
		removeContigs_if(g, LongerThanX(g, seen, opt::maxLen));

	// Remove contigs with low mean k-mer coverage.
	if (opt::minCoverage > 0)
		removeContigs_if(g, CoverageLessThan(g, seen, opt::minCoverage));

	// Remove contigs with high mean k-mer coverage.
	if (opt::maxCoverage > 0)
		removeContigs_if(g, std::not1(CoverageLessThan(g, seen, opt::maxCoverage)));

	// Remove inconsistent edges of spaceseeds
	if (argc - optind == 1) {
		const char* contigsPath(argv[optind++]);
		Contigs& contigs = g_contigs;
		if (opt::verbose > 0)
			cerr << "Reading `" << contigsPath << "'...\n";
		FastaReader in(contigsPath, FastaReader::NO_FOLD_CASE);
		for (FastaRecord rec; in >> rec;) {
			if (g_contigNames.count(rec.id) == 0)
				continue;
			assert(contigs.size() == get(g_contigNames, rec.id));
			contigs.push_back(rec.seq);
		}
		assert(in.eof());

		removeEdges_if(g, is_edge_inconsistent(g));
	}

	if (opt::verbose > 0) {
		cerr << "Graph stats after:\n";
		printGraphStats(cerr, g);
	}

	// Output the updated adjacency graph.
	if (!opt::graphPath.empty()) {
		ofstream fout(opt::graphPath.c_str());
		assert_good(fout, opt::graphPath);
		write_graph(fout, g, PROGRAM, commandLine);
		assert_good(fout, opt::graphPath);
	}

	// Assemble unambiguous paths. These need to be assembled by
	// MergeContigs before being processed by other applications.
	if (opt::assemble) {
		size_t numContigs = num_vertices(g) / 2;
		typedef vector<ContigPath> ContigPaths;
		ContigPaths paths;
		if (opt::ss)
			assemble_stranded(g, back_inserter(paths));
		else
			assemble(g, back_inserter(paths));
		g_contigNames.unlock();
		for (ContigPaths::const_iterator it = paths.begin(); it != paths.end(); ++it) {
			ContigNode u(numContigs + it - paths.begin(), false);
			string name = createContigName();
			put(vertex_name, g, u, name);
			cout << name << '\t' << *it << '\n';
		}
		g_contigNames.lock();
	}

	return 0;
}
