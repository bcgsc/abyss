#include "dialign.h"
#include "config.h"
#include "Common/Options.h"
#include "FastaReader.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include "alignGlobal.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "abyss-align"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FASTA]...\n"
"Align multiple sequences globally using either Needleman-Wunsch\n"
"or DIALIGN-TX. Groups of sequences may be separated using `#.'\n"
"\n"
" Arguments:\n"
"\n"
"  FASTA  sequences in FASTA format\n"
"\n"
" Options:\n"
"\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
" DIALIGN-TX options:\n"
"\n"
"  -D, --dialign-d=N     dialign debug level, default: 0\n"
"  -M, --dialign-m=FILE  score matrix, default: dna_matrix.scr\n"
"  -P, --dialign-p=FILE  diagonal length probability distribution\n"
"                        default: dna_diag_prob_100_exp_550000\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static int dialign_debug;
	static string dialign_score;
	static string dialign_prob;
}

static const char shortopts[] = "v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ "dialign-d",   required_argument, NULL, 'D' },
	{ "dialign-m",   required_argument, NULL, 'M' },
	{ "dialign-p",   required_argument, NULL, 'P' },
	{ NULL, 0, NULL, 0 }
};

/** Align two sequences using the Needlman-Wunsch algorithm.
 */
static void alignPair(const string& seq0, const string& seq1,
		ostream& out)
{
	NWAlignment align;
	unsigned match = alignGlobal(seq0, seq1, align);
	float identity = (float)match / align.size();
	out << align << identity << "\n\n";
}

/** Align multiple sequences using DIALIGN-TX. */
static void alignMulti(const vector<string>& seq, ostream& out)
{
	unsigned match;
	string alignment;
	string consensus = dialign(seq, alignment, match);
	float identity = (float)match / consensus.size();
	out << alignment << consensus << '\n' << identity << "\n\n";
}

/** Align the specified sequences. */
static void align(const vector<string>& seq, ostream& out)
{
	switch (seq.size()) {
	  case 0:
		return;
	  case 1:
		out << seq.front() << '\n' << 1 << "\n\n";
		return;
	  case 2:
		return alignPair(seq[0], seq[1], out);
	  default:
		return alignMulti(seq, out);
	}
}


/** Align multiple sequences. */
static void alignFile(const char* path) {
	if (opt::verbose > 0)
		cerr << "Aligning `" << path << "'\n";
	FastaReader in(path, FastaReader::NO_FOLD_CASE);
	for (vector<string> seq; in;) {
		seq.clear();
		FastaRecord fa;
		if (in >> fa)
			seq.push_back(fa.seq);
		while (in.peek() == '>' && in >> fa)
			seq.push_back(fa.seq);
		align(seq, cout);
	}
	assert(in.eof());
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
			shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case '?': die = true; break;
		case 'D': arg >> opt::dialign_debug; break;
		case 'M': arg >> opt::dialign_score; break;
		case 'P': arg >> opt::dialign_prob; break;
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

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	// Initialize dialign.
	init_parameters();
	set_parameters_dna();
	para->DEBUG = opt::dialign_debug;
	para->SCR_MATRIX_FILE_NAME = (char*)opt::dialign_score.c_str();
	para->DIAG_PROB_FILE_NAME = (char*)opt::dialign_prob.c_str();
	initDialign();

	if (optind < argc)
		for_each(&argv[optind], &argv[argc], alignFile);
	else
		alignFile("-");

	return 0;
}
