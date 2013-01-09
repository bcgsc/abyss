/**
 * Print a kmer file. A kmer file is a serialized Google sparsehash.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */

#include "Sequence.h"
#include "SequenceCollection.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <getopt.h>
#include <iostream>
#include <sstream>

using namespace std;

#define PROGRAM "kmerprint"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2013 Shaun Jackman\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... GRAPH...\n"
"Convert a binary de Bruijn graph to plain text.\n"
"  GRAPH  the binary de Bruijn graph\n"
"\n"
"  -k, --kmer=N          length of a k-mer\n"
"      --ray             output in Ray Cloud Browser format\n"
"      --tsv             output in TSV format [default]\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static bool strands;

	/** The output field separator character. */
	static int ofsInt = '\t';
	static char ofs;
}

static const char shortopts[] = "k:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer", required_argument, NULL, 'k' },
	{ "ray", no_argument, &opt::ofsInt, ';' },
	{ "tsv", no_argument, &opt::ofsInt, '\t' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static void print(const SequenceCollectionHash::value_type& seq)
{
	const KmerData& data = seq.second;
	cout << seq.first.str()
		<< opt::ofs << data.getMultiplicity()
		<< opt::ofs << data.getExtension(SENSE)
		<< opt::ofs << data.getExtension(ANTISENSE)
		<< '\n';
}

static void print(const SequenceCollectionHash::value_type& seq,
		extDirection sense)
{
	const KmerData& data = seq.second;
	cout << (sense ? reverseComplement(seq.first).str()
			: seq.first.str())
		<< opt::ofs << data.getMultiplicity(sense)
		<< opt::ofs << data.getExtension(sense)
		<< '\n';
}

static void printFile(const char* path)
{
	SequenceCollectionHash c;
	c.load(path);
	for (SequenceCollectionHash::const_iterator it = c.begin();
			it != c.end(); ++it) {
		if (it->second.deleted())
			continue;
		if (opt::strands) {
			print(*it, SENSE);
			print(*it, ANTISENSE);
		} else
			print(*it);
	}
}

int main(int argc, char* argv[])
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		  case '?':
			  exit(EXIT_FAILURE);
		  case 'k':
			  arg >> opt::k;
			  break;
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
	opt::ofs = opt::ofsInt;

	if (opt::k <= 0) {
		cerr << PROGRAM ": " << "missing -k,--kmer option\n";
		die = true;
	}

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	Kmer::setLength(opt::k);

	for_each(argv + optind, argv + argc, printFile);
	return 0;
}
