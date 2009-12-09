/** Written by Shaun Jackman <sjackman@bcgsc.ca>. */

#include "config.h"
#include "Common/Options.h"
#include "DataLayer/Options.h"
#include "Kmer.h"
#include <algorithm>
#include <climits> // for INT_MAX
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

using namespace std;

#define PROGRAM "ABYSS"

namespace opt {

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... FILE...\n"
"Assemble the input files, FILE, which may be in FASTA, FASTQ,\n"
"qseq or export format and compressed with gz, bz2 or xz.\n"
"\n"
"      --chastity                 discard unchaste reads [default]\n"
"                                 for qseq- and export-formatted files only\n"
"      --no-chastity              do not discard unchaste reads\n"
"      --trim-masked              trim masked bases from the ends of reads\n"
"                                 [default]\n"
"      --no-trim-masked           do not trim masked bases from the ends of reads\n"
"  -q, --trim-quality=THRESHOLD   trim bases from the ends of reads whose quality\n"
"                                 is less than the threshold\n"
"      --standard-fastq           zero quality is `!' (33) [default]\n"
"      --illumina-fastq           zero quality is `@' (64)\n"
"  -o, --out=FILE                 write the contigs to FILE\n"
"                                 the default is contigs.fa\n"
"  -k, --kmer=KMER_SIZE           k-mer size\n"
"  -l, --read-length=READ_LENGTH  read length\n"
"  -t, --trim-length=TRIM_LENGTH  maximum length of dangling edges to trim\n"
"  -c, --coverage=COVERAGE        remove contigs with mean k-mer coverage\n"
"                                 less than this threshold\n"
"  -b, --bubbles=N                maximum number of bubble-popping rounds\n"
"  -e, --erode=COVERAGE           erode bases at the ends of blunt contigs with\n"
"                                 coverage less than this threshold\n"
"  -E, --erode-strand=COVERAGE    erode bases at the ends of blunt contigs with\n"
"                                 coverage less than this threshold on either\n"
"                                 strand. default=1\n"
"      --coverage-hist=FILE       record the k-mer coverage histogram in FILE\n"
"  -g, --graph=FILE               generate a graph in dot format\n"
"  -s, --snp=FILE                 record popped bubbles in FILE\n"
"  -v, --verbose                  display verbose output\n"
"      --help     display this help and exit\n"
"      --version  output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

/** k-mer length */
int kmerSize = -1;

/** erosion coverage */
unsigned erode = (unsigned)-1;

/** erosion strand coverage */
unsigned erodeStrand = 1;

/** trim length */
int trimLen = -1;

/** Coverage cutoff. */
float coverage = -1;

/** Maximum number of bubble-popping rounds. */
int bubbles = INT_MAX;

/** coverage histogram path */
string coverageHistPath;

/** output contigs path */
string contigsPath = "contigs.fa";

/** temporary output contigs path
 * Each node stores its contigs in its own file temporarily.
 */
string contigsTempPath;

/** graph output */
string graphPath;

/** output bubble path */
string snpPath;

/** input FASTA files */
vector<string> inFiles;

static const char shortopts[] = "b:c:e:E:g:k:l:o:q:s:t:v";

enum { OPT_HELP = 1, OPT_VERSION, COVERAGE_HIST };

static const struct option longopts[] = {
	{ "out",         required_argument, NULL, 'o' },
	{ "kmer",        required_argument, NULL, 'k' },
	{ "read-length", required_argument, NULL, 'l' },
	{ "trim-length", required_argument, NULL, 't' },
	{ "chastity",    no_argument,       &opt::chastityFilter, 1 },
	{ "no-chastity", no_argument,       &opt::chastityFilter, 0 },
	{ "trim-masked",    no_argument,    &opt::trimMasked, 1 },
	{ "no-trim-masked", no_argument,    &opt::trimMasked, 0 },
	{ "trim-quality",   required_argument, NULL, 'q' },
	{ "standard-fastq", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-fastq", no_argument, &opt::qualityOffset, 64 },
	{ "coverage",    required_argument, NULL, 'c' },
	{ "coverage-hist", required_argument, NULL, COVERAGE_HIST },
	{ "bubbles",     required_argument, NULL, 'b' },
	{ "erode",       required_argument, NULL, 'e' },
	{ "erode-strand", required_argument, NULL, 'E' },
	{ "no-erode",    no_argument,       (int*)&erode, 0 },
	{ "graph",       required_argument, NULL, 'g' },
	{ "snp",         required_argument, NULL, 's' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Parse the specified command line. */
void parse(int argc, char* const* argv)
{
	ostringstream sargv;
	if (opt::rank <= 0) {
		char* const* last = argv + argc - 1;
		copy(argv, last, ostream_iterator<const char *>(sargv, " "));
		sargv << *last;
	}

	/** read length */
	int readLen = -1;

	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg;
		if (optarg != NULL)
			arg.str(optarg);
		switch (c) {
			case '?':
				die = true;
				break;
			case 'b':
				arg >> bubbles;
				break;
			case 'c':
				arg >> coverage;
				break;
			case 'k':
				arg >> kmerSize;
				break;
			case 'l':
				arg >> readLen;
				break;
			case COVERAGE_HIST:
				coverageHistPath = optarg;
				break;
			case 'o':
				contigsPath = optarg;
				break;
			case 'e':
				arg >> erode;
				break;
			case 'E':
				arg >> erodeStrand;
				break;
			case 't':
				arg >> trimLen;
				break;
			case 'g':
				graphPath = optarg;
				break;
			case 'q':
				arg >> opt::qualityThreshold;
				break;
			case 's':
				snpPath = optarg;
				break;
			case 'v':
				verbose++;
				break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (kmerSize <= 0) {
		cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}
	if (argv[optind] == NULL) {
		cerr << PROGRAM ": missing input sequence file argument\n";
		die = true;
	}
	if (die) {
		cerr << "Try `" PROGRAM " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	if (readLen > 0) {
		if (kmerSize > readLen) {
			cerr << PROGRAM ": k-mer size must not be larger than "
				"the read length\n";
			exit(EXIT_FAILURE);
		}
		if (trimLen < 0)
			trimLen = 6 * (readLen - kmerSize + 1);
	}

	assert(opt::qualityOffset >= 33);
	assert(opt::qualityOffset < 100);
	assert(opt::qualityThreshold <= 40);

	if (trimLen < 0)
		trimLen = kmerSize;

	Kmer::setLength(kmerSize);

	inFiles.resize(argc - optind);
	copy(&argv[optind], &argv[argc], inFiles.begin());

	if (rank >= 0) {
		ostringstream s;
		s << "contigs-" << opt::rank << ".fa";
		contigsTempPath = s.str();
	}

	if (opt::rank <= 0)
		cout << PACKAGE_STRING "\n" << sargv.str() << endl;
}

} // namespace opt
