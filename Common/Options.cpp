/** Written by Shaun Jackman <sjackman@bcgsc.ca>. */

#include "config.h"
#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

using namespace std;

namespace opt {

static const char *VERSION_MESSAGE =
PACKAGE " (ABySS) " VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2008 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PACKAGE " [OPTION]... FILE...\n"
"Assemble all input files, FILE, which may be in FASTA (.fa) format or\n"
"FASTQ format (.fastq). The contigs are written to contigs.fa.\n"
"\n"
"  -k, --kmer=KMER_SIZE           k-mer size\n"
"  -l, --read-length=READ_LENGTH  read length\n"
"  -t, --trim-length=TRIM_LENGTH  maximum length of dangling edges to trim\n"
"  -b, --bubbles=N                maximum number of bubble-popping rounds\n"
"  -e, --erode=N                  erode N bases from the ends of contigs\n"
"                                 Yields slightly shorter contigs, but fewer\n"
"                                 errors at the ends of contigs.\n"
"  -g, --graph=FILE               generate a graph in dot format\n"
"  -s, --snp=FILE                 record SNPs in FILE\n"
"  -v, --verbose                  display verbose output\n"
"      --help     display this help and exit\n"
"      --version  output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

/** k-mer length */
int kmerSize = -1;

/** read length */
int readLen = -1;

/** erode length */
int erode = -1;

/** trim length */
int trimLen = -1;

/** Maximum number of bubble-popping rounds. */
int bubbles = 5;

/** graph output */
std::string graphPath;

/** output SNP path */
std::string snpPath;

/** output SNP file */
FILE* snpFile;

/** verbose output */
int verbose = 0;

/** input FASTA files */
vector<std::string> inFiles;

static const char *shortopts = "b:e:g:k:l:s:t:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "read-length", required_argument, NULL, 'l' },
	{ "trim-length", required_argument, NULL, 't' },
	{ "bubbles",     required_argument, NULL, 'b' },
	{ "erode",       required_argument, NULL, 'e' },
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
	bool die = false;
	char c;
	while ((c = getopt_long(argc, argv, shortopts, longopts, NULL))
			!= -1) {
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
			case 'k':
				arg >> kmerSize;
				break;
			case 'l':
				arg >> readLen;
				break;
			case 'e':
				arg >> erode;
				break;	
			case 't':
				arg >> trimLen;
				break;
			case 'g':
				graphPath = optarg;
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
		cerr << PACKAGE ": " << "missing -k,--kmer option\n";
		die = true;
	}
	if (readLen <= 0) {
		cerr << PACKAGE ": " << "missing -l,--read-length option\n";
		die = true;
	}
	if (argv[optind] == NULL) {
		cerr << PACKAGE ": " << "missing input sequence file argument\n";
		die = true;
	}
	if (die) {
		cerr << "Try `" << PACKAGE
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
	if (kmerSize > readLen) {
		cerr << PACKAGE ": "
			<< "k-mer size must be smaller than the read length\n";
		exit(EXIT_FAILURE);
	}

	inFiles.resize(argc - optind);
	copy(&argv[optind], &argv[argc], inFiles.begin());

	unsigned n = readLen - kmerSize + 1;
	if (erode < 0)
		erode = n;
	if (trimLen < 0)
		trimLen = 6 * n;

	if (snpPath.length() > 0) {
		snpFile = fopen(snpPath.c_str(), "w");
		if (snpFile == NULL) {
			perror(snpPath.c_str());
			exit(EXIT_FAILURE);
		}
	}
}

} // namespace opt
