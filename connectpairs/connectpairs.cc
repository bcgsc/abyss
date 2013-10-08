/**
 * Connect pairs using a Bloom filter de Bruijn graph
 */

#include "DBGBloom.h"
#include "DBGBloomAlgorithms.h"

#include "Common/Options.h"
#include "DataLayer/FastaInterleave.h"
#include "DataLayer/Options.h"
#include "Graph/DotIO.h"
#include "Graph/Options.h"

#include <cassert>
#include <getopt.h>
#include <iostream>

#undef USESEQAN

#if USESEQAN
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/align_split.h>
#endif

using namespace std;
#if USESEQAN
using namespace seqan;
#endif

#define PROGRAM "abyss-connectpairs"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman, Hamid Mohamadi, Anthony Raymond and\n"
"Ben Vandervalk.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [READS1 READS2]...\n"
"Connect the pairs READS1 and READS2 and close the gap using\n"
"a Bloom filter de Bruijn graph.\n"
"\n"
" Options:\n"
"\n"
"  -k, --kmer=N            the size of a k-mer\n"
"      --chastity          discard unchaste reads [default]\n"
"      --no-chastity       do not discard unchaste reads\n"
"      --trim-masked       trim masked bases from the ends of reads\n"
"      --no-trim-masked    do not trim masked bases from the ends\n"
"                          of reads [default]\n"
"  -q, --trim-quality=N    trim bases from the ends of reads whose\n"
"                          quality is less than the threshold\n"
"      --standard-quality  zero quality is `!' (33)\n"
"                          default for FASTQ and SAM files\n"
"      --illumina-quality  zero quality is `@' (64)\n"
"                          default for qseq and export files\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** The size of a k-mer. */
	unsigned k;
}

/** Counters */
struct {
	size_t noPath;
	size_t uniquePath;
	size_t multiplePaths;
	size_t tooManyPaths;
} g_count;

static const char shortopts[] = "k:q:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",             required_argument, NULL, 'k' },
	{ "chastity",         no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",      no_argument, &opt::chastityFilter, 0 },
	{ "trim-masked",      no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",   no_argument, &opt::trimMasked, 0 },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "verbose",          no_argument, NULL, 'v' },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "version",          no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Load the bloom filter. */
template <typename It>
void loadBloomFilter(DBGBloom& g, It first, It last)
{
	if (opt::verbose > 0)
		cerr << "Constructing the Bloom filter\n";
	for (It it = first; it != last; ++it) {
		std::string path = *it;
		if (opt::verbose > 0)
			cerr << "Reading `" << path << "'...\n";
		g.open(path);
	}
	if (opt::verbose > 0)
		cerr << "Loaded " << num_vertices(g) << " k-mer\n";
}

#if USESEQAN
const string r1 =
"AGAATCAACCAACCGTTCAATGATATAATCAAGAGCGATATTGTAATCTTTGTTTCT";
const string r2 =
"CGACGTCCACCAATTCGTCCCTGTGCACGAGCAGTTTCCAGTCCAGCTTTTGTTCGT";
const string ins =
"AGAATCAACCAACCGTTCAATGATATAATCAAGAGCGATATTGTAATCTTTGTTTCTGTCACCCGGCCCCCACGACTCAAGGATTAGACCATAAACACCATCCTCTTCACCTATCGAACACTCAGCTTTCAGTTCAATTCCATTATTATCAAAAACATGCATAATATTAATCTTTAATCAATTTTTCACGACAATACTACTTTTATTGATAAAATTGCAACAAGTTGCTGTTGTTTTACTTTCTTTTGTACACAAAGTGTCTTTAACTTTATTTATCCCCTGCAGGAAACCTCTTATACAAAGTTGACACACCAACATCATAGATAATCGCCACCTTCTGGCGAGGAGTTCCTGCTGCAATTAATCGTCCAGCTTGTGCCCATTGTTCTGGTGTAAGTTTGGGACGACGTCCACCAATTCGTCCCTGTGCACGAGCAGTTTCCAGTCCAGCTTTTGTTCGT";

static void seqanTests()
{
	typedef String<Dna> DS;
	typedef Align<DS> Alignment;

    //DS seq1 = "TTGT";
    //DS seq2 = "TTAGT";
	DS ref = ins;
	DS seq1 = r1;
	DS seq2 = r2;

    Alignment align1;
	resize(rows(align1), 2);
	assignSource(row(align1, 0), ref);
	assignSource(row(align1, 1), seq1);
    Alignment align2;
	resize(rows(align2), 2);
	assignSource(row(align2, 0), ref);
	assignSource(row(align2, 1), seq2);

	Score<int> scoring(2, -2, -50, -100);

	cout << splitAlignment(align1, align2, scoring) << endl;
	cout << align1 << endl;
	cout << align2 << endl;

	cout << localAlignment(align1, scoring) << endl;
	cout << align1 << endl;

	cout << localAlignment(align2, scoring) << endl;
	cout << align2 << endl;
}
#endif

/** Connect a read pair. */
static void connectPair(const DBGBloom& g,
	const FastaRecord& read1, const FastaRecord& read2)
{
	const unsigned maxNumPaths = 2;
	const unsigned maxPathLen = 1000;

	vector<FastaRecord> paths;
	PathSearchResult result
		= connectPairs(read1, read2, g, paths,
				maxNumPaths, maxPathLen);
	switch (result) {
	  case NO_PATH:
		assert(paths.empty());
		++g_count.noPath;
		break;
	  case FOUND_PATH:
		assert(!paths.empty());
		if (paths.size() == 1) {
			++g_count.uniquePath;
			cout << paths.front();
		} else
			++g_count.multiplePaths;
		break;
	  case TOO_MANY_PATHS:
		++g_count.tooManyPaths;
		break;
	}
}

/** Connect read pairs. */
static void connectPairs(const DBGBloom& g, FastaInterleave& in)
{
	for (FastaRecord a, b; in >> a >> b;)
		connectPair(g, a, b);
}

/**
 * Connect pairs using a Bloom filter de Bruijn graph
 */
int main(int argc, char** argv)
{
	bool die = false;

	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			die = true; break;
		  case 'k':
			arg >> opt::k; break;
		  case 'q':
			arg >> opt::qualityThreshold; break;
		  case 'v':
			opt::verbose++; break;
		  case OPT_HELP:
			cerr << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		  case OPT_VERSION:
			cerr << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (opt::k == 0) {
		cerr << PROGRAM ": missing mandatory option `-k'\n";
		die = true;
	}

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	Kmer::setLength(opt::k);

#if USESEQAN
	seqanTests();
#endif

	DBGBloom g(opt::k);
	loadBloomFilter(g, argv + optind, argv + argc);

	if (opt::verbose > 0)
		cerr << "Connecting read pairs\n";
	FastaInterleave in(argv + optind, argv + argc,
			FastaReader::FOLD_CASE);
	connectPairs(g, in);
	assert(in.eof());

	if (opt::verbose > 0) {
		size_t n = g_count.uniquePath + g_count.noPath
			+ g_count.multiplePaths + g_count.tooManyPaths;
		cerr <<
			"Total number of read pairs: " << n << "\n"
			"No path: " << g_count.noPath
				<< " (" << setprecision(3) << (float)100
					* g_count.noPath / n << "%)\n"
			"Unique path: " << g_count.uniquePath
				<< " (" << setprecision(3) << (float)100
					* g_count.uniquePath / n << "%)\n"
			"Multiple paths: " << g_count.multiplePaths
				<< " (" << setprecision(3) << (float)100
					* g_count.multiplePaths / n << "%)\n"
			"Too many paths: " << g_count.tooManyPaths
				<< " (" << setprecision(3) << (float)100
					* g_count.tooManyPaths / n << "%)\n";
	}

	cout.flush();
	assert_good(cout, "stdout");

	return 0;
}
