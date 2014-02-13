#include "Common/Options.h"
#include "DataLayer/Options.h"
#include "ContigNode.h"
#include "ContigProperties.h"
#include "FastaReader.h"
#include "Iterator.h"
#include "Kmer.h"
#include "StringUtil.h"
#include "SuffixArray.h"
#include "Uncompress.h"
#include "UnorderedMap.h"
#include "Graph/ContigGraph.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

#define PROGRAM "AdjList"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Find overlaps of [m,k) bases. Contigs may be read from FILE(s)\n"
"or standard input. Output is written to standard output.\n"
"Overlaps of exactly k-1 bases are found using a hash table.\n"
"Overlaps of fewer than k-1 bases are found using a suffix array.\n"
"\n"
" Options:\n"
"\n"
"  -k, --kmer=K          find overlaps of up to K-1 bases\n"
"  -m, --min-overlap=M   require a minimum overlap of M bases [50]\n"
"      --adj             output the results in adj format [default]\n"
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

	/** The minimum required amount of overlap. */
	static unsigned minOverlap = 50;
}

static const char shortopts[] = "k:m:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",    required_argument, NULL, 'k' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "adj",     no_argument,       &opt::format, ADJ },
	{ "dot",     no_argument,       &opt::format, DOT },
	{ "sam",     no_argument,       &opt::format, SAM },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** A contig adjacency graph. */
typedef DirectedGraph<ContigProperties, Distance> DG;
typedef ContigGraph<DG> Graph;

/** Parse and return the coverage from the specified FASTA comment. */
static unsigned getCoverage(const string& comment)
{
	istringstream ss(comment);
	unsigned length, coverage = 0;
	ss >> length >> coverage;
	return coverage;
}

/** Add the overlaps of vseq to the graph. */
static void addOverlapsSA(Graph& g, const SuffixArray& sa,
		ContigNode v, const string& vseq)
{
	assert(!vseq.empty());
	set<ContigNode> seen;
	typedef SuffixArray::const_iterator It;
	for (string q(vseq, 0, vseq.size() - 1);
			q.size() >= opt::minOverlap; chop(q)) {
		pair<It, It> range = sa.equal_range(q);
		for (It it = range.first; it != range.second; ++it) {
			ContigNode u(it->second);
			if (seen.insert(u).second) {
				// Add the longest overlap between two vertices.
				unsigned overlap = it->first.size();
				add_edge(u, v, -overlap, static_cast<DG&>(g));
			}
		}
	}
}

/** Add overlaps of fewer than k-1 bp to the graph. */
static void addOverlapsSA(Graph& g, const vector<Kmer>& prefixes)
{
	// Construct a suffix array of the blunt contigs.
	typedef pair<string, ContigNode> Suffix;
	typedef deque<Suffix> Suffixes;
	Suffixes suffixes;
	SuffixArray sa(opt::minOverlap);

	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef graph_traits<Graph>::vertex_iterator Vit;
	pair<Vit, Vit> vertices = g.vertices();
	for (Vit it = vertices.first; it != vertices.second; ++it) {
		ContigNode u(*it);
		if (out_degree(u, g) > 0)
			continue;
		size_t uci = get(vertex_index, g,
				get(vertex_complement, g, u));
		assert(uci < prefixes.size());
		string suffix(reverseComplement(prefixes[uci]).str());
		suffixes.push_back(Suffix(suffix, u));
		sa.insert(suffixes.back());
	}
	sa.construct();

	for (Suffixes::const_iterator it = suffixes.begin();
			it != suffixes.end(); ++it) {
		V uc = get(vertex_complement, g, it->second);
		addOverlapsSA(g, sa, uc, reverseComplement(it->first));
	}
}

/** An index of suffixes of k-1 bp. */
#if defined(HAVE_STD_HASH) || defined(HAVE_STD_TR1_HASH)
typedef unordered_map<Kmer, vector<ContigNode>, hash<Kmer> >
	SuffixMap;
#else
typedef unordered_map<Kmer, vector<ContigNode>, boost::hash<Kmer> >
	SuffixMap;
#endif

/** Read contigs. Add contig properties to the graph. Add prefixes to
 * the collection and add suffixes to their index.
 */
static void readContigs(const string& path,
		Graph& g, vector<Kmer>& prefixes, SuffixMap& suffixMap)
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

		// Add the prefix to the collection.
		unsigned overlap = opt::k - 1;
		assert(seq.length() > overlap);
		Kmer prefix(seq.substr(0, overlap));
		Kmer suffix(seq.substr(seq.length() - overlap));
		prefixes.push_back(prefix);
		prefixes.push_back(reverseComplement(suffix));

		// Add the suffix to the index.
		ContigProperties vp(seq.length(), getCoverage(rec.comment));
		ContigNode u = add_vertex(vp, g);
		put(vertex_name, g, u, rec.id);
		suffixMap[suffix].push_back(u);
		suffixMap[reverseComplement(prefix)].push_back(
				get(vertex_complement, g, u));
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
			case 'm': arg >> opt::minOverlap; break;
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

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	if (opt::minOverlap == 0)
		opt::minOverlap = opt::k - 1;
	opt::minOverlap = min(opt::minOverlap, opt::k - 1);

	opt::trimMasked = false;
	Kmer::setLength(opt::k - 1);
	Graph g;
	vector<Kmer> prefixes;
	SuffixMap suffixMap(prefixes.size());
	if (optind < argc) {
		for (; optind < argc; optind++)
			readContigs(argv[optind], g, prefixes, suffixMap);
	} else
		readContigs("-", g, prefixes, suffixMap);
	g_contigNames.lock();

	// Add the overlap edges of exactly k-1 bp.
	typedef graph_traits<Graph>::vertex_descriptor V;
	if (opt::verbose > 0)
		cerr << "Finding overlaps of exactly k-1 bp...\n";
	for (vector<Kmer>::const_iterator it = prefixes.begin();
			it != prefixes.end(); ++it) {
		ContigNode v(it - prefixes.begin());
		const SuffixMap::mapped_type& edges = suffixMap[*it];
		for (SuffixMap::mapped_type::const_iterator
				itu = edges.begin(); itu != edges.end(); ++itu) {
			V uc = get(vertex_complement, g, *itu);
			V vc = get(vertex_complement, g, v);
			add_edge(vc, uc, -(int)opt::k + 1, static_cast<DG&>(g));
		}
	}
	SuffixMap().swap(suffixMap);

	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	if (opt::minOverlap < opt::k - 1) {
		// Add the overlap edges of fewer than k-1 bp.
		if (opt::verbose > 0)
			cerr << "Finding overlaps of fewer than k-1 bp...\n";
		addOverlapsSA(g, prefixes);
		if (opt::verbose > 0)
			printGraphStats(cerr, g);
	}

	// Output the graph.
	write_graph(cout, g, PROGRAM, commandLine);
	assert(cout.good());

	return 0;
}
