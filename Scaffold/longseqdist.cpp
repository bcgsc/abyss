#include "config.h"
#include "IOUtil.h"
#include "ContigNode.h"
#include "Uncompress.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphAlgorithms.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include "ContigProperties.h"
#include "SAM.h"
#include <iostream>
#include <sstream>
#include <string>
#include <getopt.h>
#include <boost/tuple/tuple.hpp>

using namespace std;
using namespace boost;

#define PROGRAM "abyss-longseqdist"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Anthony Raymond.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k<kmer> [OPTION]... SAM >DIST\n"
"Generate distance estimates between all contigs a single\n"
"read maps to.\n"
"\n"
" Arguments:\n"
"\n"
"  SAM       BWA-MEM alignments of long sequences to the assembly\n"
"  DIST      estimates of the distance between contigs\n"
"\n"
" Options:\n"
"\n"
"  -k, --kmer=N          length of a k-mer\n"
"      --min-gap=N       minimum scaffold gap length to output [200]\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties

	/** Minimum scaffold gap length to output. */
	static int minGap = 200;

	/** Verbose output. */
	int verbose; // used by PopBubbles

	/** Output format */
	int format = DOT; // used by DistanceEst
}

static const char shortopts[] = "k:n:o:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MIN_GAP };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "min-gap",     required_argument, NULL, OPT_MIN_GAP },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** A distance estimate graph. */
typedef DirectedGraph<Length, DistanceEst> DG;
typedef ContigGraph<DG> Graph;

static void processQuery(vector<Alignment>& recs, Graph& g)
{
	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef graph_traits<Graph>::edge_descriptor E;
	typedef edge_property<Graph>::type EP;
	if (recs.size() <= 1)
		return;

	sort(recs.begin(), recs.end());
	for (vector<Alignment>::const_iterator itx = recs.begin();
			itx != recs.end(); itx++) {
		for (vector<Alignment>::const_iterator ity = itx + 1;
				ity != recs.end(); ity++) {
			// if aligned to the same contig, don't draw self edge
			if (itx->contig == ity->contig)
				continue;

			// if y is subsumed in x don't draw edge
			unsigned xqend = itx->read_start_pos + itx->align_length;
			unsigned yqend = ity->read_start_pos + ity->align_length;
			if (xqend >= yqend)
				continue;

			V u = find_vertex(itx->contig, itx->isRC, g);
			V v = find_vertex(ity->contig, ity->isRC, g);
			E e;
			EP ep(opt::minGap, 1, opt::minGap);
			bool found;
			tie(e, found) = edge(u, v, g);
			if (found) {
				EP& ep = g[e];
				ep.numPairs++;
			} else
				add_edge(u, v, ep, g);
		}
	}
}

static void readAlignments(istream& in, Graph& g)
{
	SAMRecord rec, prev;
	vector<Alignment> recs;
	int i = 0;
	while (prev.isUnmapped() || prev.mapq == 0)
		in >> prev;
	recs.push_back(prev);
	while (in >> rec) {
		if (rec.isUnmapped() || rec.mapq == 0)
			continue;
		if (opt::verbose > 0 && ++i % 100000 == 0)
			cerr << "Processed " << i << " good alignments...\n";
		if (rec.qname != prev.qname) {
			processQuery(recs, g);
			recs.clear();
			prev = rec;
			recs.push_back(prev);
		} else
			recs.push_back(rec);
	}
	if (opt::verbose > 0)
		cerr << "Processed " << i << " good alignments.\n";
	processQuery(recs, g);
}

int main(int argc, char** argv)
{
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
		  case 'v':
			opt::verbose++;
			break;
		  case OPT_MIN_GAP:
			arg >> opt::minGap;
			break;
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
		cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}

	if (argc - optind != 1) {
		cerr << PROGRAM ": incorrect number of arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	if (opt::verbose > 0)
		cerr << "Reading graph file '" << argv[optind] << "`...\n";
	Graph g;
	ifstream in(argv[optind]);
	in >> g;
	g_contigNames.lock();
	if (opt::verbose > 0)
		cerr << "Finished reading graph.\n";
	if (in.eof()) {
		cerr << PROGRAM ": there are no alignments\n";
		exit(EXIT_FAILURE);
	}

	if (opt::verbose > 0)
		cerr << "Processing alignments from '" << argv[optind] << "`...\n";
	readAlignments(in, g);
	write_dot(cout, g);
	return 0;
}
