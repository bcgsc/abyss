#include "Common/Options.h"
#include "DataLayer/Options.h"
#include "ContigGraph.h"
#include "ContigNode.h"
#include "ContigProperties.h"
#include "DirectedGraph.h"
#include "FastaReader.h"
#include "GraphIO.h"
#include "GraphUtil.h"
#include "HashMap.h"
#include "Iterator.h"
#include "Kmer.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

using namespace std;

#define PROGRAM "AdjList"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Find all contigs that overlap by exactly k-1 bases. Contigs may be read\n"
"from FILE(s) or standard input. Output is written to standard output.\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"      --adj             output the results in adj format [DEFAULT]\n"
"      --dot             output the results in dot format\n"
"      --sam             output the results in SAM format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by GraphIO
	int format; // used by GraphIO
}

static const char shortopts[] = "k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",    required_argument, NULL, 'k' },
	{ "adj",     no_argument,       &opt::format, ADJ },
	{ "dot",     no_argument,       &opt::format, DOT },
	{ "sam",     no_argument,       &opt::format, SAM },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** A contig adjacency graph. */
typedef DirectedGraph<ContigProperties> DG;
typedef ContigGraph<DG> Graph;

/** Return the distance between two vertices. */
static inline int get(edge_distance_t, const Graph&,
		graph_traits<Graph>::edge_descriptor)
{
	return -opt::k + 1;
}

/** The two terminal Kmer of a contig and its length and coverage. */
struct ContigEndSeq {
	unsigned length;
	unsigned coverage;
	Kmer l;
	Kmer r;
	ContigEndSeq(unsigned length, unsigned coverage,
			const Kmer& l, const Kmer& r)
		: length(length), coverage(coverage), l(l), r(r) { }
};

static unsigned getCoverage(const string& comment)
{
	istringstream ss(comment);
	unsigned length, coverage = 0;
	ss >> length >> coverage;
	return coverage;
}

static void readContigs(string path, vector<ContigEndSeq>* pContigs)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << path << "'...\n";

	unsigned count = 0;
	FastaReader in(path.c_str(), FastaReader::FOLD_CASE);
	for (FastaRecord rec; in >> rec;) {
		const Sequence& seq = rec.seq;
		if (count++ == 0) {
			// Detect colour-space contigs.
			opt::colourSpace = isdigit(seq[0]);
		} else {
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));
		}

		ContigID id(rec.id);
		assert(id == pContigs->size());

		unsigned overlap = opt::k - 1;
		assert(seq.length() > overlap);
		Kmer seql(seq.substr(seq.length() - overlap, overlap));
		Kmer seqr(seq.substr(0, overlap));
		pContigs->push_back(ContigEndSeq(seq.length(),
					getCoverage(rec.comment),
					seql, seqr));
	}
	assert(in.eof());
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
			case 'k': arg >> opt::k; break;
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

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	opt::trimMasked = false;
	Kmer::setLength(opt::k - 1);
	vector<ContigEndSeq> contigs;
	if (optind < argc) {
		for_each(argv + optind, argv + argc,
				bind2nd(ptr_fun(readContigs), &contigs));
	} else
		readContigs("-", &contigs);
	ContigID::lock();

	if (opt::verbose > 0)
		cerr << "Read " << contigs.size() << " contigs\n";

	typedef hash_map<Kmer, vector<ContigNode>, hashKmer> KmerMap;
	vector<KmerMap> ends(2, KmerMap(contigs.size()));
	for (vector<ContigEndSeq>::const_iterator it = contigs.begin();
			it != contigs.end(); ++it) {
		unsigned i = it - contigs.begin();
		ends[0][it->l].push_back(ContigNode(i, false));
		ends[1][reverseComplement(it->l)].push_back(
				ContigNode(i, true));
		ends[1][it->r].push_back(ContigNode(i, false));
		ends[0][reverseComplement(it->r)].push_back(
				ContigNode(i, true));
	}

	// Add the vertices.
	Graph g;
	for (vector<ContigEndSeq>::const_iterator it = contigs.begin();
			it != contigs.end(); ++it)
		add_vertex(ContigProperties(it->length, it->coverage), g);
	
	// Add the edges.
	typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator itu = vit.first; itu != vit.second; ++itu) {
		ContigNode u(*itu);
		const ContigEndSeq& contig = contigs[ContigID(u)];
		const Kmer& kmer = !u.sense() ? contig.l : contig.r;
		const KmerMap::mapped_type& edges = ends[!u.sense()][kmer];
		for (KmerMap::mapped_type::const_iterator
				itv = edges.begin(); itv != edges.end(); ++itv)
			add_edge<DG>(u, *itv ^ u.sense(), g);
	}

	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	write_graph(cout, g, PROGRAM, commandLine);
	assert(cout.good());

	return 0;
}
