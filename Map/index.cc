#include "config.h"
#include "FMIndex.h"
#include "IOUtil.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

using namespace std;

#define PROGRAM "abyss-index"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... FILE\n"
"Build an FM-index of FILE and store it in FILE.fm.\n"
"\n"
"      --decompress        decompress the INDEX\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** Decompress the index. */
	static int decompress;

	/** Verbose output. */
	static int verbose;
}

static const char shortopts[] = "v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "decompress", no_argument, &opt::decompress, 1 },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

int main(int argc, char **argv)
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
	}

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (argc - optind > 1) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	if (opt::decompress) {
		// Decompress the index.
		string fmPath(argv[optind]);
		if (fmPath.size() < 4
				|| !equal(fmPath.end() - 3, fmPath.end(), ".fm"))
			fmPath.append(".fm");
		string faPath(fmPath, 0, fmPath.size() - 3);

		FMIndex fmIndex(fmPath);
		ofstream out(faPath.c_str());
		fmIndex.decompress(ostream_iterator<uint8_t>(out, ""));
		assert_good(out, faPath);
		return 0;
	}

	const char* faPath(argv[optind]);
	ostringstream ss;
	ss << faPath << ".fm";
	string fmPath = ss.str();

	FMIndex f;
	f.setAlphabet("\1\nACGT");
	f.buildFmIndex(faPath, 4);

	ofstream out(fmPath.c_str());
	assert_good(out, fmPath);
	f.save(out);
	assert_good(out, fmPath);

	return 0;
}
