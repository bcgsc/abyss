#include "dialign.h"
#include "config.h"
#include "Common/Options.h"
#include "FastaReader.h"
#include "IOUtil.h"
#include "needleman_wunsch.h"
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
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FASTA]\n"
"Align multiple sequences globally using either Needleman-Wunsch\n"
"or DIALIGN-TX. Groups of sequences may be separated using `#.'\n"
"  FASTA  sequences in FASTA format\n"
"\n"
" Options:\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
" DIALIGN-TX options:\n"
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
	if (opt::verbose > 1)
		cerr << align << '\n';
	out << identity << '\n';
}

/** Align multiple sequences using DIALIGN-TX. */
static void alignMulti(const vector<string>& seq, ostream& out)
{
	unsigned match;
	string consensus = dialign(seq, match);
	float identity = (float)match / consensus.size();
	if (opt::verbose > 1)
		cerr << consensus << "\n\n";
	out << identity << '\n';
}

/** Align the specified sequences. */
static void align(const vector<string>& seq, ostream& out)
{
	switch (seq.size()) {
	  case 0:
		return;
	  case 1:
		if (opt::verbose > 1)
			cerr << seq.front() << "\n\n";
		out << 1 << '\n';
		return;
	  case 2:
		return alignPair(seq[0], seq[1], out);
	  default:
		return alignMulti(seq, out);
	}
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
	}

	if (argc - optind < 0) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	} else if (argc - optind > 1) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	const char* faPath = optind < argc ? argv[optind++] : "-";

	// Initialize dialign.
	init_parameters();
	set_parameters_dna();
	para->DEBUG = opt::dialign_debug;
	para->SCR_MATRIX_FILE_NAME = (char*)opt::dialign_score.c_str();
	para->DIAG_PROB_FILE_NAME = (char*)opt::dialign_prob.c_str();
	initDialign();

	// Align input sequences.
	FastaReader in(faPath, FastaReader::NO_FOLD_CASE);
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
	return 0;
}
