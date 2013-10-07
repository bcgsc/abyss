/**
 * Connect pairs using a Bloom filter de Bruijn graph
 */

#include "DBGBloom.h"
#include "DBGBloomAlgorithms.h"

#include "Common/Options.h"
#include "DataLayer/Options.h"
#include "Graph/DotIO.h"
#include "Graph/Options.h"
#include "Graph/GraphUtil.h"

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
"Usage: " PROGRAM " [OPTION]... READS1 READS2\n"
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
static void loadBloomFilter(DBGBloom& g, const string& path)
{
	g.open(path);
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

static void findPaths() { assert(false); abort(); }

static void mergePathts() { assert(false); abort(); }

static void buildRead() { assert(false); abort(); }

static void processRead()
{
	findPaths();
	mergePathts();
	buildRead();
}

static void processReads()
{
	// for each read
	processRead();
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

	if (argc - optind > 2) {
		cerr << PROGRAM ": too many arguments\n";
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
	assert(optind < argc);
	loadBloomFilter(g, argv[optind]);
	if (opt::verbose > 0)
		printGraphStats(cerr, g);

	write_dot(cout, g);

	processReads();

	return 0;
}
