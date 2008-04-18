#include "config.h"
#include <getopt.h>
#include <iostream>
#include <sstream>

using namespace std;

namespace opt {

unsigned kmerSize;
unsigned readLen;
string fastaFile;

static const char *shortopts = "k:l:";

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "read-length", required_argument, NULL, 'l' },
};

void parse(int argc, char* const* argv)
{
	char c;
	while ((c = getopt_long(argc, argv, shortopts, longopts, NULL))
			!= -1) {
		istringstream arg(optarg);
		switch (c) {
			case 'k':
				arg >> kmerSize;
				break;
			case 'l':
				arg >> readLen;
				break;
		}
	}

	bool die = false;
	if (kmerSize == 0) {
		cerr << PACKAGE ": " << "missing -k,--kmer option\n";
		die = true;
	}
	if (readLen == 0) {
		cerr << PACKAGE ": " << "missing -l,--read-length option\n";
		die = true;
	}
	if (argc - optind > 1) {
		cerr << PACKAGE ": " << "unexpected argument: "
			<< argv[optind+1] << "\n";
		die = true;
	}
	if (argv[optind] == NULL) {
		cerr << PACKAGE ": " << "missing FASTA file argument\n";
		die = true;
	}
	if (die)
		exit(EXIT_FAILURE);

	fastaFile = argv[optind++];
}

} // namespace opt
