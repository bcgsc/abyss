/**
 * Find contigs that overlap and end due to a lack of coverage.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "config.h"
#include "Common/Options.h"
#include "ContigGraph.h"
#include "ContigGraphAlgorithms.h"
#include "ContigProperties.h"
#include "Estimate.h"
#include "FastaReader.h"
#include "MapGraph.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "Overlap"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... CONTIGS ADJ DIST\n"
"Find overlaps between blunt contigs that have negative distance estimates.\n"
"Output the small contigs that fill in the gaps.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -m, --min=OVERLAP     require a minimum of OVERLAP bases\n"
"                        default is 5 bases\n"
"      --scaffold        join contigs with Ns [default]\n"
"      --no-scaffold     do not scaffold\n"
"      --mask-repeat     join contigs at a simple repeat and mask the repeat\n"
"                        [default]\n"
"      --no-merge-repeat don't join contigs at a repeat\n"
"  -g, --graph=FILE      write the contig adjacency graph to FILE\n"
"  -o, --out=FILE        write result to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	int k; // used by ContigGraph
	static unsigned minimum_overlap = 5;
	static int mask = 1;
	static int scaffold = 1;

	/** Write the contig adjacency graph to this file. */
	static string graphPath;

	/** Write the new contigs to this file. */
	static string out;

	int dot; // used by Estimate
}

static const char shortopts[] = "g:k:m:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",    required_argument, NULL, 'k' },
	{ "min",     required_argument, NULL, 'm' },
	{ "scaffold", no_argument,      &opt::scaffold, 1 },
	{ "no-scaffold", no_argument,   &opt::scaffold, 0 },
	{ "mask-repeat",    no_argument, &opt::mask, 1 },
	{ "no-merge-repeat", no_argument, &opt::mask, 0 },
	{ "out",     required_argument, NULL, 'o' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Contig sequences. */
static vector<string> g_contigs;

/** Contig adjacency graph. */
typedef ContigGraph<DirectedGraph<ContigProperties> > Graph;
static Graph g_graph;

static struct {
	unsigned overlap;
	unsigned scaffold;
	unsigned none;
	unsigned tooshort;
	unsigned homopolymer;
	unsigned motif;
	unsigned ambiguous;
} stats;

/** Return the sequence of the specified contig. */
static string sequence(const ContigNode& id)
{
	const string& seq = g_contigs[id.id()];
	return id.sense() ? reverseComplement(seq) : seq;
}

static unsigned findOverlap(const ContigNode& t_id,
		const ContigNode& h_id,
		bool& mask)
{
	mask = false;
	string t = sequence(t_id);
	string h = sequence(h_id);
	unsigned len = min(t.length(), h.length());
	vector<unsigned> overlaps;
	overlaps.reserve(len);
	for (unsigned overlap = len; overlap >= 1; overlap--) {
		if (t.substr(t.length()-overlap, overlap)
				== h.substr(0, overlap))
			overlaps.push_back(overlap);
	}

	if (opt::verbose > 0) {
		cout << t_id << '\t' << h_id;
		for (vector<unsigned>::const_iterator i = overlaps.begin();
				i != overlaps.end(); ++i)
			cout << '\t' << *i;
		cout << '\n';
	}

	if (overlaps.empty()) {
		stats.none++;
		return 0;
	}

	if (overlaps[0] < opt::minimum_overlap) {
		stats.tooshort++;
		return 0;
	}

	if (overlaps.size() >= 3
			&& overlaps[0]-overlaps[1] == overlaps[1]-overlaps[2]) {
		// Homopolymer run or motif.
		if (overlaps[0]-overlaps[1] == 1)
			stats.homopolymer++;
		else
			stats.motif++;
		mask = true;
	}

	return overlaps[0];
}

static FastaRecord newContig(const ContigNode& t, const ContigNode& h,
		int dist, const string& seq)
{
	ostringstream comment;
	comment << seq.length() << " 0 " << t << ' ' << h << ' ' << dist;
	return FastaRecord(ContigID::create().str(), comment.str(), seq);
}

static FastaRecord overlapContigs(const ContigNode& t_id,
		const ContigNode& h_id, unsigned overlap, bool mask)
{
	string t = sequence(t_id);
	string h = sequence(h_id);
	assert(overlap < (unsigned)opt::k - 1);
	unsigned gap = opt::k - 1 - overlap;
	string a(t, t.length() - opt::k+1, gap);
	string o(h, 0, overlap);
	string b(h, overlap, gap);
	if (mask)
		transform(o.begin(), o.end(), o.begin(), ptr_fun(::tolower));
	return newContig(t_id, h_id, -overlap, a + o + b);
}

static FastaRecord mergeContigs(
		const ContigNode& t, const ContigNode& h,
		const Estimate& est, unsigned overlap, bool mask)
{
	if (overlap > 0) {
		stats.overlap++;
		int dist = -overlap;
		int diff = dist - est.distance;
		if ((unsigned)abs(diff) > allowedError(est.stdDev))
			cerr << "warning: overlap does not agree with estimate\n"
				<< '\t' << t << ',' << h << ' '
				<< dist << ' ' << est << endl;
		return overlapContigs(t, h, overlap, mask);
	} else {
		assert(opt::scaffold);
		stats.scaffold++;
		if (opt::verbose > 0)
			cout << t << '\t' << h << "\t(" << est.distance << ")\n";
		string gap = est.distance <= 0 ? string("n")
			: string(est.distance, 'N');
		const string& ts = sequence(t);
		const string& hs = sequence(h);
		unsigned overlap = opt::k - 1;
		return newContig(t, h, est.distance,
				ts.substr(ts.length() - overlap) + gap
				+ hs.substr(0, overlap));
	}
}

struct Overlap {
	const Estimate est;
	const unsigned overlap;
	const bool mask;
	Overlap(const Estimate& est, unsigned overlap, bool mask)
		: est(est), overlap(overlap), mask(mask) { }
};

static FastaRecord mergeContigs(
		const ContigNode& t, const ContigNode& h,
		const Overlap& overlap)
{
	return mergeContigs(t, h, overlap.est, overlap.overlap,
			overlap.mask);
}

/** The scaffold graph. Edges join two blunt contigs that are joined
 * by a distance estimate. */
typedef map<ContigNode, set<ContigNode> > OverlapGraph;
static OverlapGraph g_scaffoldGraph;

/** A subgraph of g_scaffoldGraph containing only the edges that
 * overlap (are not scaffolded). */
static OverlapGraph g_overlapGraph;

/** The amount of overlap (attributes of OverlapGraph). */
typedef pair<ContigNode, ContigNode> Edge;
typedef map<Edge, Overlap> OverlapGraphAttr;
static OverlapGraphAttr g_overlaps;

static void removeVertex(OverlapGraph& g, const ContigNode& u)
{
	typedef graph_traits<OverlapGraph>::adjacency_iterator
		adjacency_iterator;
	pair<adjacency_iterator, adjacency_iterator>
		adj = adjacent_vertices(u, g);
	for (adjacency_iterator it = adj.first; it != adj.second; ++it)
		remove_edge(~*it, ~u, g);
	clear_vertex(u, g);
	remove_vertex(u, g);
}

static void findOverlap(
		ContigID refID, bool rc, const Estimate& est)
{
	if (refID == est.contig.id()
			|| (est.distance >= 0 && !opt::scaffold))
		return;
	ContigNode ref(refID, false);
	const ContigNode& pair = est.contig;
	const ContigNode& t = rc ? pair : ref;
	const ContigNode& h = rc ? ref : pair;
	if (g_graph.out_degree(t) > 0 || g_graph.in_degree(h) > 0)
		return;

	bool mask = false;
	unsigned overlap
		= est.distance - (int)allowedError(est.stdDev) <= 0
		? findOverlap(t, h, mask) : 0;
	if (mask && !opt::mask)
		return;
	if (overlap > 0) {
		g_overlapGraph[t].insert(h);
		g_overlapGraph[~h].insert(~t);
	}
	if (overlap > 0 || opt::scaffold) {
		g_scaffoldGraph[t].insert(h);
		g_scaffoldGraph[~h].insert(~t);
		// Store the edge attribute only once.
		if (g_overlaps.count(Edge(~h, ~t)) == 0)
			g_overlaps.insert(make_pair(Edge(t, h),
						Overlap(est, overlap, mask)));
	}
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

static void readContigs(const char *contigPath)
{
	FastaReader in(contigPath, FastaReader::FOLD_CASE);
	for (FastaRecord rec; in >> rec;)
		g_contigs.push_back(rec.seq);
	assert(in.eof());
	assert(!g_contigs.empty());
	opt::colourSpace = isdigit(g_contigs[0][0]);
}

int main(int argc, char *const argv[])
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'g': arg >> opt::graphPath; break;
			case 'k': arg >> opt::k; break;
			case 'm': arg >> opt::minimum_overlap; break;
			case 'o': arg >> opt::out; break;
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

	if (opt::out.empty()) {
		cerr << PROGRAM ": " << "missing -o,--out option\n";
		die = true;
	}

	if (argc - optind < 3) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (argc - optind > 3) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	const char* contigPath(argv[optind++]);
	string adjPath(argv[optind++]);
	string estPath(argv[optind++]);

	readContigs(contigPath);

	// Read the contig adjacency graph.
	ifstream fin(adjPath.c_str());
	assert_open(fin, adjPath);
	fin >> g_graph;
	assert(fin.eof());

	ofstream out(opt::out.c_str());
	assert(out.is_open());
	ifstream in(estPath.c_str());
	assert_open(in, estPath);

	for (EstimateRecord er; in >> er;) {
		for (int rc = false; rc <= true; ++rc) {
			const vector<Estimate>& ests = er.estimates[rc];
			for (EstimateVector::const_iterator iter = ests.begin();
					iter != ests.end(); ++iter)
				findOverlap(er.refID, rc, *iter);
		}
	}
	assert(in.eof());
	in.close();
	ContigID::unlock();

	if (opt::verbose > 1)
		cout << "digraph overlap {\n"
			<< dot_writer<OverlapGraph>(g_scaffoldGraph)
			<< "}\n";

	// First, give priority to overlapping edges (not scaffolded).
	for (OverlapGraphAttr::const_iterator it = g_overlaps.begin();
			it != g_overlaps.end(); ++it) {
		const ContigNode& t = source(it->first, g_overlapGraph),
			  h = target(it->first, g_overlapGraph);
		if (it->second.overlap == 0) {
			// This edge is scaffolded.
		} else if (contiguous_out(g_overlapGraph, t)) {
			assert(*g_overlapGraph[t].begin() == h);
			FastaRecord contig = mergeContigs(t, h, it->second);
			out << contig;
			assert(out.good());

			// Add the new contig to the adjacency graph.
			Graph::vertex_descriptor v = g_graph.add_vertex(
					ContigProperties(contig.seq.length(), 0));
			g_graph.add_edge(t, v);
			g_graph.add_edge(v, h);

			// Remove the vertices incident to this edge from the
			// scaffold graph.
			removeVertex(g_scaffoldGraph, t);
			removeVertex(g_scaffoldGraph, ~h);
		} else
			stats.ambiguous++;
	}

	// Second, handle scaffolded edges.
	for (OverlapGraphAttr::const_iterator it = g_overlaps.begin();
			it != g_overlaps.end(); ++it) {
		const ContigNode& t = source(it->first, g_scaffoldGraph),
			  h = target(it->first, g_scaffoldGraph);
		if (it->second.overlap > 0) {
			// This edge is not scaffolded.
		} else if (g_scaffoldGraph[t].count(h) == 0) {
			// This edge involved a vertex that has already been used
			// and removed.
		} else if (contiguous_out(g_scaffoldGraph, t)) {
			assert(*g_scaffoldGraph[t].begin() == h);
			FastaRecord contig = mergeContigs(t, h, it->second);
			out << contig;
			assert(out.good());

			// Add the new contig to the adjacency graph.
			Graph::vertex_descriptor v = g_graph.add_vertex(
					ContigProperties(contig.seq.length(), 0));
			g_graph.add_edge(t, v);
			g_graph.add_edge(v, h);
		} else
			stats.ambiguous++;
	}
	out.close();

	if (!opt::graphPath.empty()) {
		// Output the updated adjacency graph.
		ofstream fout(opt::graphPath.c_str());
		assert(fout.good());
		fout << g_graph;
		assert(fout.good());
	}

	cout << "Overlap: " << stats.overlap << "\n"
		"Scaffold: " << stats.scaffold << "\n"
		"No overlap: " << stats.none << "\n"
		"Insignificant (<" << opt::minimum_overlap << "bp): "
		<< stats.tooshort << "\n"
		"Homopolymer: " << stats.homopolymer << "\n"
		"Motif: " << stats.motif << "\n"
		"Ambiguous: " << stats.ambiguous << "\n";
	return 0;
}
