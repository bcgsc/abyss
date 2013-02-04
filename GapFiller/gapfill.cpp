#include "SAM.h"
#include "FastaReader.h"

#include "config.h"
#include <cstdlib>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


using namespace std;

#define PROGRAM "abyss-gapfill"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Anthony Raymond.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

//TODO
static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Write read pairs that map to the same contig to the file SAME.\n"
"Write read pairs that map to different contigs to stdout.\n"
"Alignments may be in FILE(s) or standard input.\n"
"\n"
"      --no-qname        set the qname to * [default]\n"
"      --qname           do not alter the qname\n"
"  -l, --min-align=N     the minimal alignment size [1]\n"
"  -s, --same=SAME       write properly-paired reads to this file\n"
"  -h, --hist=FILE       write the fragment size histogram to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

//TODO
namespace opt {
	static string alignPath;
	static int verbose;
}

//TODO
static const char shortopts[] = "h:l:s:v";

//TODO
enum { OPT_HELP = 1, OPT_VERSION };

//TODO
static const struct option longopts[] = {
	{ "min-align", required_argument, NULL, 'l' },
	{ "hist",    required_argument, NULL, 'h' },
	{ "same",    required_argument, NULL, 's' },
	{ "verbose", no_argument,       NULL, 'v' },
	{ "help",    no_argument,       NULL, OPT_HELP },
	{ "version", no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

typedef map<string, FastaRecord> Scaffolds

static void readContigs(string path, Scaffolds& scaffs)
{
	FastaReader in(path.c_str());
	FastaRecord rec;
	while (in >> rec) {
		size_t i = rec.seq.find_first_of("N")
		if (i != string::npos) {
			
		}
	}
}


int main(int argc, char* const* argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
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

	if (argc - optind != 1) {
		cerr << PROGRAM ": incorrect number of arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	readScaffolds();
	readAlignments();

	fillGaps();
	writeContigs();

	return 0;
}
