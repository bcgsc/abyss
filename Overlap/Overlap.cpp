/**
 * Find contigs that overlap and end due to a lack of coverage.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "config.h"
#include "Common/Options.h"
#include "ContigProperties.h"
#include "Estimate.h"
#include "FastaReader.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/ref.hpp>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits> // for UINT_MAX
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace boost::lambda;

#define PROGRAM "Overlap"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k<kmer> -o<out.fa> [OPTION]... CONTIGS ADJ DIST\n"
"Find overlaps between blunt contigs that have negative distance\n"
"estimates. Add edges to the overlap graph.\n"
"\n"
" Options:\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -m, --min=OVERLAP     require a minimum of OVERLAP bases\n"
"                        default is 5 bases\n"
"      --scaffold        join contigs with Ns [default]\n"
"      --no-scaffold     do not scaffold\n"
"      --mask-repeat     join contigs at a simple repeat and mask\n"
"                        the repeat sequence [default]\n"
"      --no-merge-repeat don't join contigs at a repeat\n"
"      --SS              expect contigs to be oriented correctly\n"
"      --no-SS           no assumption about contig orientation [default]\n"
"  -g, --graph=FILE      write the contig adjacency graph to FILE\n"
"  -o, --out=FILE        write result to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigGraph
	static unsigned minimum_overlap = 5;
	static int mask = 1;
	static int scaffold = 1;

	/** Run a strand-specific RNA-Seq assembly. */
	static int ss;

	/** The acceptable error of a distance estimate. */
	unsigned distanceError = 6;

	/** Write the contig adjacency graph to this file. */
	static string graphPath;

	/** Write the new contigs to this file. */
	static string out;

 	/** Output format */
 	int format = ADJ; // used by ContigProperties
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
	{ "SS",            no_argument,       &opt::ss, 1 },
	{ "no-SS",         no_argument,       &opt::ss, 0 },
	{ "out",     required_argument, NULL, 'o' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Contig sequences. */
static vector<string> g_contigs;

/** Contig adjacency graph. */
typedef ContigGraph<DirectedGraph<ContigProperties, Distance> > Graph;

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

static unsigned findOverlap(const Graph& g,
		const ContigNode& t_id,
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
		cout << get(vertex_name, g, t_id)
			<< '\t' << get(vertex_name, g, h_id);
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

static FastaRecord newContig(const Graph& g,
		const ContigNode& t, const ContigNode& v,
		int dist, const string& seq)
{
	ostringstream comment;
	comment << seq.length() << " 0 "
		<< get(vertex_name, g, t) << ' '
		<< get(vertex_name, g, v) << ' ' << dist;
	return FastaRecord(createContigName(), comment.str(), seq);
}

/** An overlap of two sequences. */
struct Overlap : public DistanceEst {
	unsigned overlap;
	bool mask;

	Overlap() : overlap(UINT_MAX), mask(false) { }
	Overlap(int) { assert(false); }
	Overlap(const DistanceEst& est, unsigned overlap, bool mask)
		: DistanceEst(est), overlap(overlap), mask(mask) { }

	bool operator==(const Overlap& o) const
	{
		return overlap == o.overlap;
	}

	operator Distance() const
	{
		assert(overlap > 0);
		return -overlap;
	}

	friend ostream& operator<<(ostream& out, const Overlap& o)
	{
		return out << "d="
			<< (o.overlap > 0 ? -(int)o.overlap : o.distance);
	}
};

/** Create a contig representing the gap between contigs u and v. */
static FastaRecord createGapContig(const Graph& g,
		const ContigNode& u, const ContigNode& v,
		const Overlap& o)
{
	assert(opt::scaffold);
	assert(o.overlap == 0);
	stats.scaffold++;
	int distance = o.distance;
	if (opt::verbose > 0)
		cout << get(vertex_name, g, u)
			<< '\t' << get(vertex_name, g, v)
			<< "\t(" << distance << ")\n";
	assert(distance < 100000);
	string gap = distance <= 0 ? string("n")
		: string(distance, 'N');
	const string& useq = sequence(u);
	const string& vseq = sequence(v);
	unsigned overlap = opt::k - 1; // by convention
	return newContig(g, u, v, distance,
			useq.substr(useq.length() - overlap) + gap
			+ vseq.substr(0, overlap));
}

/** The scaffold graph. Edges join two blunt contigs that are joined
 * by a distance estimate. */
typedef ContigGraph<DirectedGraph<NoProperty, Overlap> > OverlapGraph;

/**
 * Check for an overlap between the specified pair of contigs.
 * Add the size of the overlap to the edge properties. Add the
 * complementary edge if it does not exist in the graph.
 * @param goverlap the contig overlap graph
 * @param g the scaffold graph
 * @return true if the contigs overlap
 */
static bool checkEdgeForOverlap(const Graph& goverlap,
		OverlapGraph& g,
		graph_traits<OverlapGraph>::edge_descriptor e)
{
	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef edge_bundle_type<OverlapGraph>::type EP;

	V u = source(e, g), v = target(e, g);
	V uc = get(vertex_complement, g, u);
	V vc = get(vertex_complement, g, v);
	assert(u != v);
	assert(u != vc);
	EP& ep = g[e];
	if (ep.overlap != UINT_MAX) {
		// Found the complementary overlap.
		return ep.overlap > 0 || opt::scaffold;
	}
	if (ep.distance >= 0 && !opt::scaffold) {
		// Positive distance estimate and not scaffolding.
		return false;
	}
	if (out_degree(u, goverlap) > 0 || in_degree(v, goverlap) > 0) {
		// Not blunt.
		return false;
	}

	bool mask = false;
	unsigned overlap
		= ep.distance - (int)allowedError(ep.stdDev) <= 0
		? findOverlap(goverlap, u, v, mask) : 0;
	if (mask && !opt::mask) {
		// Ambiguous overlap.
		return false;
	}
	if (overlap == 0 && !opt::scaffold) {
		// No overlap and not scaffolding.
		return false;
	}
	ep.overlap = overlap;
	ep.mask = mask;
	pair<E, bool> ecomplement = edge(vc, uc, g);
	if (ecomplement.second) {
		// Modify the complementary edge.
		g[ecomplement.first] = ep;
	} else {
		// Add the complementary edge.
		assert(vc != u);
		add_edge(vc, uc, ep,
				static_cast<OverlapGraph::base_type&>(g));
	}
	return true;
}

static void findOverlap(const Graph& g,
		ContigID refID, bool rc,
		const ContigNode& pair,
		const DistanceEst& est,
		OverlapGraph& out)
{
	if (refID == pair.id()
			|| (est.distance >= 0 && !opt::scaffold))
		return;
	ContigNode ref(refID, false);
	const ContigNode& t = rc ? pair : ref;
	const ContigNode& h = rc ? ref : pair;
	if (out_degree(t, g) > 0 || in_degree(h, g) > 0
			|| edge(t, h, out).second)
		return;

	bool mask = false;
	unsigned overlap
		= est.distance - (int)allowedError(est.stdDev) <= 0
		? findOverlap(g, t, h, mask) : 0;
	if (mask && !opt::mask)
		return;
	if (overlap > 0 || opt::scaffold)
		add_edge(t, h, Overlap(est, overlap, mask), out);
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
	assert_good(fin, adjPath);
	Graph graph;
	fin >> graph;
	assert(fin.eof());
	g_contigNames.lock();

	// Open the output file.
	ofstream out(opt::out.c_str());
	assert_good(out, opt::out);

	// Read the scaffold graph.
	ifstream in(estPath.c_str());
	assert_good(in, estPath);

	// Find overlapping contigs.
	OverlapGraph scaffoldGraph(graph.num_vertices() / 2);
	if (in.peek() == 'd') {
		// dot graph format
		in >> scaffoldGraph;
		assert(in.eof());
		if (opt::verbose > 0)
			printGraphStats(cout, scaffoldGraph);
		remove_edge_if(
				!boost::lambda::bind(checkEdgeForOverlap,
					boost::cref(graph), boost::ref(scaffoldGraph),
					_1),
				static_cast<OverlapGraph::base_type&>(scaffoldGraph));
	} else {
		// dist graph format
		for (EstimateRecord er; in >> er;) {
			for (int sense = false; sense <= true; ++sense) {
				typedef vector<
					pair<ContigNode, DistanceEst> > Estimates;
				const Estimates& ests = er.estimates[sense];
				for (Estimates::const_iterator it = ests.begin();
						it != ests.end(); ++it)
					findOverlap(graph, er.refID, sense,
							it->first, it->second,
							scaffoldGraph);
			}
		}
		assert(in.eof());
	}
	in.close();

	if (opt::verbose > 1)
		cout << dot_writer(scaffoldGraph);

	typedef graph_traits<OverlapGraph>::vertex_descriptor
		vertex_descriptor;
	typedef graph_traits<OverlapGraph>::vertex_iterator
		vertex_iterator;
	typedef graph_traits<OverlapGraph>::edge_descriptor
		edge_descriptor;
	typedef graph_traits<OverlapGraph>::out_edge_iterator
		out_edge_iterator;

	/** The overlapping edges (d<0) of scaffoldGraph. */
	OverlapGraph overlapGraph(num_vertices(graph) / 2);

	/** The canonical edges of scaffoldGraph. */
	unsigned numOverlaps = num_edges(scaffoldGraph) / 2;
	typedef vector<edge_descriptor> Edges;
	Edges edges;
	edges.reserve(numOverlaps);

	// Create the set of canonical edges and the overlap subgraph.
	std::pair<vertex_iterator, vertex_iterator>
		uit = vertices(scaffoldGraph);
	for (vertex_iterator u = uit.first; u != uit.second; ++u) {
		std::pair<out_edge_iterator, out_edge_iterator>
			vit = out_edges(*u, scaffoldGraph);
		for (out_edge_iterator e = vit.first; e != vit.second; ++e) {
			vertex_descriptor v = target(*e, scaffoldGraph);
			assert(*u != v);
			if (v < *u)
				continue;
			edges.push_back(*e);
			const Overlap& ep = get(edge_bundle, scaffoldGraph, e);
			if (ep.overlap > 0)
				add_edge(*u, v, ep, overlapGraph);
		}
	}
	assert(edges.size() == numOverlaps);

	// First, give priority to overlapping edges (not scaffolded).
	for (Edges::const_iterator it = edges.begin();
			it != edges.end(); ++it) {
		const ContigNode& t = source(*it, overlapGraph),
			  h = target(*it, overlapGraph);
		if (!edge(t, h, overlapGraph).second) {
			// This edge is scaffolded.
			continue;
		}
		const Overlap& overlap = get(edge_bundle, overlapGraph, *it);
		assert(overlap.overlap > 0);
		if (contiguous_out(overlapGraph, t)) {
			stats.overlap++;
			assert(*adjacent_vertices(t, overlapGraph).first == h);
			add_edge(t, h, overlap, graph);

			// Clear the out-edges of t and the in-edges of h.
			clear_out_edges(t, scaffoldGraph);
			clear_in_edges(h, scaffoldGraph);
		} else
			stats.ambiguous++;
	}
	overlapGraph.clear();

	// Second, handle scaffolded edges.
	g_contigNames.unlock();
	for (Edges::const_iterator it = edges.begin();
			it != edges.end(); ++it) {
		const ContigNode& t = source(*it, scaffoldGraph),
			  h = target(*it, scaffoldGraph);
		if (!edge(t, h, scaffoldGraph).second) {
			// This edge involved a vertex that has already been used
			// and removed.
			continue;
		}
		const Overlap& overlap = get(edge_bundle,
				scaffoldGraph, *it);
		if (overlap.overlap > 0) {
			// This edge is not scaffolded.
		} else if (contiguous_out(scaffoldGraph, t)) {
			assert(*adjacent_vertices(t, scaffoldGraph).first == h);
			ContigNode t1 = t, h1 = h;
			if (opt::ss && t.sense() && h.sense()) {
				t1 = h ^ true;
				h1 = t ^ true;
			}
			FastaRecord contig = createGapContig(graph,
					t1, h1, overlap);
			out << contig;
			assert(out.good());

			// Add the new contig to the adjacency graph.
			vertex_descriptor v = add_vertex(
					ContigProperties(contig.seq.length(), 0), graph);
			put(vertex_name, graph, v, contig.id);
			add_edge(t1, v, graph);
			add_edge(v, h1, graph);
		} else
			stats.ambiguous++;
	}
	g_contigNames.lock();
	out.close();

	if (!opt::graphPath.empty()) {
		// Output the updated adjacency graph.
		ofstream fout(opt::graphPath.c_str());
		assert_good(fout, opt::graphPath);
		write_graph(fout, graph, PROGRAM, commandLine);
		assert_good(fout, opt::graphPath);
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
