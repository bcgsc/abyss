/** Written by Shaun Jackman <sjackman@bcgsc.ca>. */

#include "config.h"
#include "Common/Options.h"
#include "DataLayer/Options.h"
#include <algorithm>
#include <cassert>
#include <climits> // for INT_MAX
#include <getopt.h>
#include <iostream>
#include <vector>
#include "DataBase/Options.h"
#include <sstream>
#include <iterator>

using namespace std;

#define PROGRAM "ABYSS"

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

namespace opt {

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson, Shaun Jackman and Anthony Raymond.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " -k<kmer> -o<output.fa> [OPTION]... FILE...\n"
"Assemble the input files, FILE, which may be in FASTA, FASTQ,\n"
"qseq, export, SAM or BAM format and compressed with gz, bz2 or xz.\n"
"\n"
" Options:\n"
"\n"
"      --chastity        discard unchaste reads [default]\n"
"      --no-chastity     do not discard unchaste reads\n"
"      --trim-masked     trim masked bases from the ends of reads\n"
"                        [default]\n"
"      --no-trim-masked  do not trim masked bases from the ends of\n"
"                        reads\n"
"  -q, --trim-quality=N  trim bases from the ends of reads whose\n"
"                        quality is less than the threshold\n"
"  -Q, --mask-quality=N  mask all low quality bases as `N'\n"
"  --standard-quality    zero quality is `!' (33)\n"
"                        default for FASTQ and SAM files\n"
"  --illumina-quality    zero quality is `@' (64)\n"
"                        default for qseq and export files\n"
"      --SS              assemble in strand-specific mode\n"
"      --no-SS           do not assemble in strand-specific mode\n"
"  -o, --out=FILE        write the contigs to FILE\n"
"  -k, --kmer=N          the length of a k-mer (when -K is not set) [<=" STR(MAX_KMER) "]\n"
"                        or the span of a k-mer pair (when -K is set)\n"
"  -K, --single-kmer=N   the length of a single k-mer in a k-mer pair\n"
"  -t, --trim-length=N   maximum length of blunt contigs to trim [k]\n"
"  -c, --coverage=FLOAT  remove contigs with mean k-mer coverage\n"
"                        less than this threshold\n"
"      --kc=N            remove all k-mers with multiplicity < N [0]\n"
"  -b, --bubbles=N       pop bubbles shorter than N bp [3*k]\n"
"  -b0, --no-bubbles     do not pop bubbles\n"
"  -e, --erode=N         erode bases at the ends of blunt contigs with coverage\n"
"                        less than this threshold [round(sqrt(median))]\n"
"  -E, --erode-strand=N  erode bases at the ends of blunt contigs\n"
"                        with coverage less than this threshold on\n"
"                        either strand [1 if sqrt(median) > 2 else 0]\n"
"  --coverage-hist=FILE  write the k-mer coverage histogram to FILE\n"
"  -m, --mask-cov        do not include kmers containing masked bases in\n"
"                        coverage calculations [experimental]\n"
"  -s, --snp=FILE        record popped bubbles in FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"      --db=FILE         specify path of database repository in FILE\n"
"      --library=NAME    specify library NAME for database\n"
"      --strain=NAME     specify strain NAME for database\n"
"      --species=NAME    specify species NAME for database\n"
"\n"
" ABYSS Options: (won't work with ABYSS-P)\n"
"\n"
"  -g, --graph=FILE      generate a graph in dot format\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

/** The length of a single k-mer.
 * When singleKmerSize == 0, paired dBG is disabled [default].
 * When singleKmerSize == -1, paired dBG is enabled and is an error.
 * The caller of opt::parse should set singleKmerSize = -1 to enable paired dBG.
 */
int singleKmerSize = 0;

/** The length of a k-mer pair, including the gap. */
int kmerSize = -1;
int k; // used by Graph

/** k-mer range */
int kMin = -1;
int kMax = -1;
int kStep = 1;

/** erosion coverage */
unsigned erode = (unsigned)-1;

/** erosion strand coverage */
unsigned erodeStrand = (unsigned)-1;

/** trim length */
int trimLen = -1;

/** Coverage cutoff. */
float coverage = -1;

/** Minimum k-mer multiplicity cutoff. */
unsigned kc = 0;

/** Pop bubbles shorter than N bp. */
int bubbleLen = -1;

/** Whether to run a strand-specific assembly. */
int ss = 0;

/**
 * do not include kmers containing masked bases in
 * coverage calculations (experimental)
 */
bool maskCov = false;

/** coverage histogram path */
string coverageHistPath;

/** output contigs path */
string contigsPath;

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

string db;
vector<string> metaVars;

/** commandline specific to assembly */
string assemblyCmd;

static const char shortopts[] = "b:c:e:E:g:k:K:mo:Q:q:s:t:v";

enum { OPT_HELP = 1, OPT_VERSION, COVERAGE_HIST, OPT_DB, OPT_LIBRARY, OPT_STRAIN, OPT_SPECIES, OPT_KC };

static const struct option longopts[] = {
	{ "out",         required_argument, NULL, 'o' },
	{ "kmer",        required_argument, NULL, 'k' },
	{ "single-kmer", required_argument, NULL, 'K' },
	{ "trim-length", required_argument, NULL, 't' },
	{ "chastity",    no_argument,       &opt::chastityFilter, 1 },
	{ "no-chastity", no_argument,       &opt::chastityFilter, 0 },
	{ "trim-masked",    no_argument,    &opt::trimMasked, 1 },
	{ "no-trim-masked", no_argument,    &opt::trimMasked, 0 },
	{ "trim-quality",   required_argument, NULL, 'q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "SS",          no_argument,       &opt::ss, 1 },
	{ "no-SS",       no_argument,       &opt::ss, 0 },
	{ "coverage",    required_argument, NULL, 'c' },
	{ "kc",          required_argument, NULL, OPT_KC },
	{ "coverage-hist", required_argument, NULL, COVERAGE_HIST },
	{ "bubble-length", required_argument, NULL, 'b' },
	{ "no-bubbles",  no_argument,       &opt::bubbleLen, 0 },
	{ "erode",       required_argument, NULL, 'e' },
	{ "erode-strand", required_argument, NULL, 'E' },
	{ "no-erode",    no_argument,       (int*)&erode, 0 },
	{ "mask-cov",    no_argument, NULL, 'm' },
	{ "graph",       required_argument, NULL, 'g' },
	{ "snp",         required_argument, NULL, 's' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ "db",          required_argument, NULL, OPT_DB },
	{ "library",     required_argument, NULL, OPT_LIBRARY },
	{ "strain",      required_argument, NULL, OPT_STRAIN },
	{ "species",     required_argument, NULL, OPT_SPECIES },
	{ NULL, 0, NULL, 0 }
};

int getVvalue()
{
	return opt::verbose;
}

string getUvalue()
{
	return opt::db;
}

string getCommand()
{
	return opt::assemblyCmd;
}

vector<string> getMetaValue()
{
	return opt::metaVars;
}

/** Parse the specified command line. */
void parse(int argc, char* const* argv)
{
	ostringstream sargv;
	if (opt::rank <= 0) {
		char* const* last = argv + argc - 1;
		copy(argv, last, ostream_iterator<const char *>(sargv, " "));
		sargv << *last;
	}

	opt::assemblyCmd = sargv.str();
	opt::metaVars.resize(3);
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?':
				die = true;
				break;
			case 'b':
				arg >> bubbleLen;
				break;
			case 'c':
				arg >> coverage;
				break;
			case 'k':
				arg >> kmerSize;
				k = kmerSize;
				kMin = kmerSize;
				switch (arg.get()) {
				  case ',':
					arg >> kMax;
					kStep = kMax - kMin;
					break;
				  case '-':
					arg >> kMax;
					if (arg.get() == ':')
						arg >> kStep;
					break;
				  default:
					kMax = kmerSize;
				}
				assert(kMin <= kMax);
				break;
			case 'K':
				arg >> singleKmerSize;
				break;
			case COVERAGE_HIST:
				getline(arg, coverageHistPath);
				break;
			case 'm':
				maskCov = true;
				break;
			case 'o':
				getline(arg, contigsPath);
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
				getline(arg, graphPath);
				break;
			case 'q':
				arg >> opt::qualityThreshold;
				break;
			case 'Q':
				arg >> opt::internalQThreshold;
				break;
			case 's':
				getline(arg, snpPath);
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
			case OPT_DB:
				arg >> db;
				break;
			case OPT_LIBRARY:
				arg >> opt::metaVars[0];
				break;
			case OPT_STRAIN:
				arg >> opt::metaVars[1];
				break;
			case OPT_SPECIES:
				arg >> opt::metaVars[2];
				break;
			case OPT_KC:
				arg >> opt::kc;
				break;
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	if (singleKmerSize < 0) {
		cerr << PROGRAM ": missing -K,--single-kmer option\n";
		die = true;
	}
	if (kmerSize <= 0) {
		cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}
	if (contigsPath.empty()) {
		cerr << PROGRAM ": missing -o,--out option\n";
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

	assert(opt::qualityThreshold <= 40);
	assert(opt::internalQThreshold <= 40);

	if (opt::rank <= 0
			&& opt::coverage >= 0 && opt::erode == (unsigned)-1)
		cerr << "warning: -c,--coverage was specified, "
			"but -e,--erode was not specified\n"
			"Previously, the default was -e2 (or --erode=2)." << endl;

	if (trimLen < 0)
		trimLen = kmerSize;
	if (bubbleLen < 0)
		bubbleLen = 3*kmerSize;
	assert(bubbleLen == 0 || bubbleLen > kmerSize);
	if (bubbleLen == 0)
		snpPath.clear();

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
