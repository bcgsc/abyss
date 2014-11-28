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
#include <set>
#include <vector>

#if _SQL
#include "DataBase/Options.h"
#include "DataBase/DB.h"
#endif

#if PAIRED_DBG
#include "PairedDBG/KmerPair.h"
#endif

using namespace std;

#define PROGRAM "AdjList"

#if _SQL
DB db;
#endif
static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k<kmer> [OPTION]... [FILE]...\n"
"Find overlaps of [m,k) bases. Contigs may be read from FILE(s)\n"
"or standard input. Output is written to standard output.\n"
"Overlaps of exactly k-1 bases are found using a hash table.\n"
"Overlaps of fewer than k-1 bases are found using a suffix array.\n"
"\n"
" Options:\n"
"\n"
"  -k, --kmer=N          The length of a k-mer pair, including the gap\n"
"  -K, --single-kmer=N   The length of a single k-mer\n"
"  -m, --min-overlap=M   require a minimum overlap of M bases [50]\n"
"      --adj             output the graph in ADJ format [default]\n"
"      --asqg            output the graph in ASQG format\n"
"      --dot             output the graph in GraphViz format\n"
"      --gv              output the graph in GraphViz format\n"
"      --gfa             output the graph in GFA format\n"
"      --sam             output the graph in SAM format\n"
"      --SS              expect contigs to be oriented correctly\n"
"      --no-SS           no assumption about contig orientation\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
#if _SQL
"      --db=FILE         specify path of database repository in FILE\n"
"      --library=NAME    specify library NAME for database\n"
"      --strain=NAME     specify strain NAME for database\n"
"      --species=NAME    specify species NAME for database\n"
#endif
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
#if _SQL
	string url;
	dbVars metaVars;
#endif
	unsigned k; // used by GraphIO

	/** Length of a single-kmer in a kmer pair */
	unsigned singleKmerSize = 0;

	int format; // used by GraphIO

	/** Run a strand-specific RNA-Seq assembly. */
	static int ss;

	/** The minimum required amount of overlap. */
	static unsigned minOverlap = 50;
}

static const char shortopts[] = "k:K:m:v";

#if _SQL
enum { OPT_HELP = 1, OPT_VERSION, OPT_DB, OPT_LIBRARY, OPT_STRAIN, OPT_SPECIES };
#else
enum { OPT_HELP = 1, OPT_VERSION };
#endif

static const struct option longopts[] = {
	{ "kmer",    required_argument, NULL, 'k' },
	{ "single-kmer", required_argument, NULL, 'K' },
	{ "min-overlap", required_argument, NULL, 'm' },
	{ "adj",     no_argument,       &opt::format, ADJ },
	{ "asqg",    no_argument,       &opt::format, ASQG },
	{ "dot",     no_argument,       &opt::format, DOT },
	{ "gv",      no_argument,       &opt::format, DOT },
	{ "gfa",     no_argument,       &opt::format, GFA },
	{ "sam",     no_argument,       &opt::format, SAM },
	{ "SS",      no_argument,       &opt::ss, 1 },
	{ "no-SS",   no_argument,       &opt::ss, 0 },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
#if _SQL
	{ "db",      required_argument, NULL, OPT_DB },
	{ "library", required_argument, NULL, OPT_LIBRARY },
	{ "strain",  required_argument, NULL, OPT_STRAIN },
	{ "species", required_argument, NULL, OPT_SPECIES },
#endif
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
			if (opt::ss && u.sense() != v.sense())
				continue;
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
	typedef vector<Suffix> Suffixes;
	Suffixes suffixes;

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
	}

	SuffixArray sa(opt::minOverlap);
	for (Suffixes::const_iterator it = suffixes.begin();
			it != suffixes.end(); ++it)
		sa.insert(*it);
	sa.construct();

	for (Suffixes::const_iterator it = suffixes.begin();
			it != suffixes.end(); ++it) {
		V uc = get(vertex_complement, g, it->second);
		addOverlapsSA(g, sa, uc, reverseComplement(it->first));
	}
}

/** Read contigs. Add contig properties to the graph. Add prefixes to
 * the collection and add suffixes to their index.
 */
template <class KmerType>
static void readContigs(const string& path,
		Graph& g, vector<KmerType>& prefixes,
		unordered_map<KmerType, vector<ContigNode> >& suffixMap)
{
	typedef unordered_map<KmerType, vector<ContigNode> > SuffixMap;

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
		KmerType prefix(seq.substr(0, overlap));
		KmerType suffix(seq.substr(seq.length() - overlap));
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

/** Build contig overlap graph. */
template <class KmerType>
static void buildOverlapGraph(Graph& g, vector<KmerType>& prefixes,
	unordered_map<KmerType, vector<ContigNode> >& suffixMap)
{
	// Add the overlap edges of exactly k-1 bp.

	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef unordered_map<KmerType, vector<ContigNode> > SuffixMap;
	typedef typename vector<KmerType>::const_iterator PrefixIterator;
	typedef const typename SuffixMap::mapped_type Edges;
	typedef typename SuffixMap::mapped_type::const_iterator SuffixIterator;

	if (opt::verbose > 0)
		cerr << "Finding overlaps of exactly k-1 bp...\n";
	for (PrefixIterator it = prefixes.begin(); it != prefixes.end(); ++it) {
		ContigNode v(it - prefixes.begin());
		Edges& edges = suffixMap[*it];
		for (SuffixIterator itu = edges.begin(); itu != edges.end(); ++itu) {
			V uc = get(vertex_complement, g, *itu);
			V vc = get(vertex_complement, g, v);
			if (opt::ss && uc.sense() != vc.sense())
				continue;
			add_edge(vc, uc, -(int)opt::k + 1, static_cast<DG&>(g));
		}
	}
	SuffixMap().swap(suffixMap);

	if (opt::verbose > 0)
		printGraphStats(cerr, g);
}

template <class KmerType>
void loadDataStructures(Graph& g, vector<KmerType>& prefixes,
	unordered_map<KmerType, vector<ContigNode> >& suffixMap,
	int argc, char** argv)
{
	if (optind < argc) {
		for (; optind < argc; optind++)
			readContigs(argv[optind], g, prefixes, suffixMap);
	} else
		readContigs("-", g, prefixes, suffixMap);
	g_contigNames.lock();
}

/** Build contig overlap graph for standard de Bruijn graph */
void buildOverlapGraph(Graph& g, int argc, char** argv)
{
	Kmer::setLength(opt::k - 1);

	vector<Kmer> prefixes;
	unordered_map<Kmer, vector<ContigNode> >
		suffixMap(prefixes.size());

	loadDataStructures(g, prefixes, suffixMap, argc, argv);
	buildOverlapGraph(g, prefixes, suffixMap);

	if (opt::minOverlap < opt::k - 1) {
		// Add the overlap edges of fewer than k-1 bp.
		if (opt::verbose > 0)
			cerr << "Finding overlaps of fewer than k-1 bp...\n";
		addOverlapsSA(g, prefixes);
		if (opt::verbose > 0)
			printGraphStats(cerr, g);
	}
}

/** Build contig overlap graph for paired de Bruijn graph */
void buildPairedOverlapGraph(Graph& g, int argc, char** argv)
{
	Kmer::setLength(opt::singleKmerSize - 1);
	KmerPair::setLength(opt::k - 1);

	vector<KmerPair> prefixes;
	unordered_map<KmerPair, vector<ContigNode> >
		suffixMap(prefixes.size());

	loadDataStructures(g, prefixes, suffixMap, argc, argv);
	buildOverlapGraph(g, prefixes, suffixMap);
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

#if _SQL
	opt::metaVars.resize(3);
#endif

	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'K': arg >> opt::singleKmerSize; break;
			case 'm': arg >> opt::minOverlap; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
#if _SQL
			case OPT_DB:
				arg >> opt::url; break;
			case OPT_LIBRARY:
				arg >> opt::metaVars[0]; break;
			case OPT_STRAIN:
				arg >> opt::metaVars[1]; break;
			case OPT_SPECIES:
				arg >> opt::metaVars[2]; break;
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

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	if (opt::minOverlap == 0)
		opt::minOverlap = opt::k - 1;
	opt::minOverlap = min(opt::minOverlap, opt::k - 1);

#if _SQL
	init (db, opt::url, opt::verbose, PROGRAM, opt::getCommand(argc, argv), opt::metaVars);
#endif
	opt::trimMasked = false;

	// contig overlap graph
	Graph g;

#if PAIRED_DBG
	if (opt::singleKmerSize > 0)
		buildPairedOverlapGraph(g, argc, argv);
	else
		buildOverlapGraph(g, argc, argv);
#else
	buildOverlapGraph(g, argc, argv);
#endif

	// Output the graph.
	write_graph(cout, g, PROGRAM, commandLine);
	assert(cout.good());
#if _SQL
	vector<int> vals = make_vector<int>()
		<< opt::ss
		<< opt::k;
	vector<int> new_vals = passGraphStatsVal(g);
	vals.insert(vals.end(), new_vals.begin(), new_vals.end());

	vector<string> keys = make_vector<string>()
		<< "SS"
		<< "K"
		<< "V"
		<< "E"
		<< "degree0pctg"
		<< "degree1pctg"
		<< "degree234pctg"
		<< "degree5pctg"
		<< "degree_max";

	for (unsigned i=0; i<vals.size(); i++)
		addToDb(db, keys[i], vals[i]);
#endif

	return 0;
}
