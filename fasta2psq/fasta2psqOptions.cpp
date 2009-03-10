/** Written by Jared Simpson <jsimpson@bcgsc.ca>. */

#include "config.h"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>

using namespace std;

namespace pp_opt {

static const char *VERSION_MESSAGE =
PACKAGE " (ABySS) " VERSION "\n"
"Written by Jared Simpson.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: PreProcess [OPTION]... FASTA_FILE\n"
"PreProcess FASTA_FILE.\n"
"\n"
"  -o, --outfile=<file>        output filename\n"
"      --help     display this help and exit\n"
"      --version  output version information and exit\n"
"\n"
"Report bugs to " PACKAGE_BUGREPORT "\n";

/** input FASTA path */
string outFile;
string fastaFile;

static const char *shortopts = "o:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "outfile",     required_argument, NULL, 'o' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Parse the specified command line. */
void parse(int argc, char* const* argv)
{
	char c;
	while ((c = getopt_long(argc, argv, shortopts, longopts, NULL))
			!= -1) {
		istringstream arg;
		if (optarg != NULL)
			arg.str(optarg);
		switch (c) {
			case 'o':
				arg >> outFile;
				break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	bool die = false;
	
	if (argc - optind > 1) {
		cerr << PACKAGE ": " << "unexpected argument: "
			<< argv[optind+1] << "\n";
		die = true;
	}
	if (argv[optind] == NULL) {
		cerr << PACKAGE ": " << "missing FASTA file argument\n";
		die = true;
	}
	if (outFile.empty()) {
		cerr << PACKAGE ": " << "missing -o,--outfile option\n";
		die = true;
	}	
	if (die) {
		cerr << "Try `" << PACKAGE
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}
	
	fastaFile = argv[optind++];
		
	// Make sure the correct extension is set
	const char *PACKED_SEQ_EXT = ".psq";
	if(outFile.substr(outFile.length() - 4) != PACKED_SEQ_EXT)
	{
		cerr << PACKAGE ": " << "invalid extension on output: " << PACKED_SEQ_EXT << " is required\n";
		exit(EXIT_FAILURE);
	}
}

} // namespace opt
