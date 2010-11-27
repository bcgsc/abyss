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
#include <cassert>
#include <cctype>
#include <cstdlib>
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
"Find overlaps of exactly k-1 bases. Contigs may be read from\n"
"FILE(s) or standard input. Output is written to standard output.\n"
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

static unsigned getCoverage(const string& comment)
{
	istringstream ss(comment);
	unsigned length, coverage = 0;
	ss >> length >> coverage;
	return coverage;
}

/** Read contigs. Add contig properties to the graph. Add prefixes to
 * the collection.
 */
static void readContigs(const string& path,
		Graph& g, vector<Kmer>& prefixes)
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

		ContigID::insert(rec.id);
		ContigProperties vp(seq.length(), getCoverage(rec.comment));
		add_vertex(vp, g);

		unsigned overlap = opt::k - 1;
		assert(seq.length() > overlap);
		prefixes.push_back(Kmer(seq.substr(0, overlap)));
		prefixes.push_back(reverseComplement(
					Kmer(seq.substr(seq.length() - overlap))));
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
	Graph g;
	vector<Kmer> prefixes;
	if (optind < argc) {
		for (; optind < argc; optind++)
			readContigs(argv[optind], g, prefixes);
	} else
		readContigs("-", g, prefixes);
	ContigID::lock();

	// Index the suffixes using a hash table.
	typedef hash_map<Kmer, vector<ContigNode>, hashKmer> KmerMap;
	KmerMap suffixes(KmerMap(prefixes.size()));
	for (vector<Kmer>::const_iterator it = prefixes.begin();
			it != prefixes.end(); ++it) {
		ContigNode u(it - prefixes.begin());
		suffixes[reverseComplement(*it)].push_back(~u);
	}

	// Add the edges.
	for (vector<Kmer>::const_iterator it = prefixes.begin();
			it != prefixes.end(); ++it) {
		ContigNode v(it - prefixes.begin());
		const KmerMap::mapped_type& edges = suffixes[*it];
		for (KmerMap::mapped_type::const_iterator
				itu = edges.begin(); itu != edges.end(); ++itu)
			add_edge<DG>(~v, ~*itu, g);
	}

	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	// Output the graph.
	write_graph(cout, g, PROGRAM, commandLine);
	assert(cout.good());

	return 0;
}
