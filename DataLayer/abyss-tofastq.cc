/** Convert various file formats to FASTQ format.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 */
#include "config.h"
#include "DataLayer/Options.h"
#include "FastaInterleave.h"
#include "FastaReader.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>

using namespace std;

#define PROGRAM "abyss-tofastq"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Convert to FASTQ format. The input format may be FASTA, FASTQ,\n"
"qseq, export, SAM or BAM format and compressed with gz, bz2 or xz\n"
"and may be tarred.\n"
"\n"
"      --cat               concatenate the records [default]\n"
"  -i, --interleave        interleave the records\n"
"      --fastq             ouput FASTQ format [default]\n"
"      --fasta             ouput FASTA format\n"
"      --chastity          discard unchaste reads [default]\n"
"                          for qseq, export and SAM files only\n"
"      --no-chastity       do not discard unchaste reads\n"
"      --trim-masked       trim masked bases from the ends of reads\n"
"      --no-trim-masked    do not trim masked bases from the ends\n"
"                          of reads [default]\n"
"  -q, --trim-quality=N    trim bases from the ends of reads whose\n"
"                          quality is less than the threshold\n"
"      --standard-quality  zero quality is `!' (33)\n"
"                          default for FASTQ and SAM files\n"
"      --illumina-quality  zero quality is `@' (64)\n"
"                          default for qseq and export files\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** Interleave the lines of each file. */
	static int interleave;

	static int toFASTQ = 1;
	static int verbose;
}

static const char shortopts[] = "iq:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "cat",              no_argument, &opt::interleave, 0 },
	{ "interleave",       no_argument, &opt::interleave, 1 },
	{ "fasta",            no_argument, &opt::toFASTQ, 0 },
	{ "fastq",            no_argument, &opt::toFASTQ, 1 },
	{ "chastity",         no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity",      no_argument, &opt::chastityFilter, 0 },
	{ "trim-masked",      no_argument, &opt::trimMasked, 1 },
	{ "no-trim-masked",   no_argument, &opt::trimMasked, 0 },
	{ "trim-quality",     required_argument, NULL, 'q' },
	{ "standard-quality", no_argument, &opt::qualityOffset, 33 },
	{ "illumina-quality", no_argument, &opt::qualityOffset, 64 },
	{ "help",             no_argument, NULL, OPT_HELP },
	{ "version",          no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** FastaReader flags. */
static const int FASTAREADER_FLAGS
	= FastaReader::NO_FOLD_CASE | FastaReader::CONVERT_QUALITY;

/** Total count. */
static struct {
	unsigned records;
	unsigned characters;
} g_total;

template <class Record>
static void convert(const char* path)
{
	FastaReader in(path, FASTAREADER_FLAGS);
	unsigned records = 0, characters = 0;
	for (Record record; in >> record;) {
		cout << record;
		assert_good(cout, "stdout");
		records++;
		characters += record.seq.size();
	}
	assert(in.eof());

	g_total.records += records;
	g_total.characters += characters;
	if (opt::verbose)
		cerr << records << '\t'
			<< characters << '\t'
			<< path << '\n';
}

/** Interleave the records. */
template <typename Record>
static void interleave(char** first, char** last)
{
	FastaInterleave in(first, last, FASTAREADER_FLAGS);
	unsigned records = 0, characters = 0;
	for (Record record; in >> record;) {
		cout << record;
		assert_good(cout, "stdout");
		records++;
		characters += record.seq.size();
	}
	assert(in.eof());

	g_total.records += records;
	g_total.characters += characters;
}

int main(int argc, char** argv)
{
	opt::trimMasked = false;

	if (string(argv[0]).find("tofasta") != string::npos)
		opt::toFASTQ = false;

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
			case 'i':
				opt::interleave = true;
				break;
			case 'q':
				arg >> opt::qualityThreshold;
				break;
			case 'v':
				opt::verbose++;
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
	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	typedef void (*F)(const char*);
	F convertFasta = convert<FastaRecord>;
	F convertFastq = convert<FastqRecord>;
	F f = opt::toFASTQ ? convertFastq : convertFasta;

	if (optind == argc) {
		f("-");
	} else if (!opt::interleave || argc - optind == 1) {
		// Concatenate.
		for_each(argv + optind, argv + argc, f);
	} else {
		// Interleave.
		if (opt::toFASTQ)
			interleave<FastqRecord>(argv + optind, argv + argc);
		else
			interleave<FastaRecord>(argv + optind, argv + argc);
	}

	if (opt::verbose && argc - optind > 1)
		cerr << g_total.records << '\t'
			<< g_total.characters << '\t'
			<< "total\n";

	cout.flush();
	assert_good(cout, "stdout");
	return 0;
}
