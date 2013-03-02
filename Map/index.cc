#include "config.h"
#include "BitUtil.h"
#include "FastaIndex.h"
#include "FMIndex.h"
#include "IOUtil.h"
#include "MemoryUtil.h"
#include "StringUtil.h"
#include "Uncompress.h"
#include <algorithm>
#include <cctype> // for toupper
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <string>

using namespace std;

#define PROGRAM "abyss-index"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2012 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... FILE\n"
"Build an FM-index of FILE and store it in FILE.fm.\n"
"\n"
" Options:\n"
"\n"
"      --both              build both FAI and FM indexes [default]\n"
"      --fai               build a FAI index\n"
"      --fm                build a FM index\n"
"      --fa2bwt            build the BWT directly without the SA\n"
"      --bwt2fm            build the FM index from the BWT\n"
"  -a, --alphabet=STRING   use the alphabet STRING [-ACGT]\n"
"  -s, --sample=N          sample the suffix array [16]\n"
"  -d, --decompress        decompress the index FILE\n"
"  -c, --stdout            write output to standard output\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** Sample the suffix array. */
	static unsigned sampleSA = 16;

	/** Which indexes to create. */
	enum { NONE, FAI, FM, BOTH };
	static int indexes = BOTH;

	/** Build the BWT directly without the SA. */
	static int fa2bwt;

	/** Build the FM index from the BWT. */
	static int bwt2fm;

	/** The alphabet. */
	static string alphabet = "-ACGT";

	/** Decompress the index. */
	static bool decompress;

	/** Write output to standard output. */
	static bool toStdout;

	/** Verbose output. */
	static int verbose;
}

static const char shortopts[] = "a:cds:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "both", no_argument, &opt::indexes, opt::BOTH },
	{ "fai", no_argument, &opt::indexes, opt::FAI },
	{ "fm", no_argument, &opt::indexes, opt::FM },
	{ "fa2bwt", no_argument, &opt::fa2bwt, true },
	{ "bwt2fm", no_argument, &opt::bwt2fm, true },
	{ "alphabet", optional_argument, NULL, 'a' },
	{ "decompress", no_argument, NULL, 'd' },
	{ "sample", required_argument, NULL, 's' },
	{ "stdout", no_argument, NULL, 'c' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Index the specified FASTA file. */
static void indexFasta(const string& faPath, const string& faiPath)
{
	cerr << "Reading `" << faPath << "'...\n";
	FastaIndex fai;
	fai.index(faPath);

	if (opt::verbose > 0)
		cerr << "Read " << fai.size() << " contigs.\n";

	cerr << "Writing `" << faiPath << "'...\n";
	ofstream out(faiPath.c_str());
	assert_good(out, faiPath);
	out << fai;
	out.flush();
	assert_good(out, faiPath);
}

/** Build the FM index from the BWT. */
static void buildFMIndexFromBWT(FMIndex& fm, const string& path)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << path << "'...\n";
	std::vector<FMIndex::value_type> bwt;
	readFile(path.c_str(), bwt);
	assert(bwt.size() > 1);

	if (opt::alphabet.empty()) {
		fm.setAlphabet(bwt.begin(), bwt.end());
		std::cerr << "The alphabet has "
			<< fm.alphabetSize() << " symbols.\n";
	} else
		fm.setAlphabet(opt::alphabet);

	fm.encode(bwt.begin(), bwt.end());
	fm.sampleSA(opt::sampleSA);
	fm.assignBWT(bwt.begin(), bwt.end());
}

/** Build an FM index of the specified file. */
static void buildFMIndex(FMIndex& fm, const string& path)
{
	if (opt::verbose > 0)
		std::cerr << "Reading `" << path << "'...\n";
	std::vector<FMIndex::value_type> s;
	readFile(path.c_str(), s);

	size_t MAX_SIZE = numeric_limits<FMIndex::sais_size_type>::max();
	if (s.size() > MAX_SIZE) {
		std::cerr << PROGRAM << ": `" << path << "', "
			<< toSI(s.size())
			<< "B, must be smaller than "
			<< toSI(MAX_SIZE) << "B\n";
		exit(EXIT_FAILURE);
	}

	// Set the alphabet.
	transform(s.begin(), s.end(), s.begin(), ::toupper);
	if (opt::alphabet.empty()) {
		fm.setAlphabet(s.begin(), s.end());
		std::cerr << "The alphabet has "
			<< fm.alphabetSize() << " symbols.\n";
	} else
		fm.setAlphabet(opt::alphabet);

	if (opt::fa2bwt) {
		// Build the BWT first.
		s.push_back(0);
		fm.buildBWT(s.begin(), s.end() - 1);
		fm.sampleSA(opt::sampleSA);
		fm.assignBWT(s.begin(), s.end());
	} else {
		// Construct the suffix array first.
		fm.assign(s.begin(), s.end());
		fm.sampleSA(opt::sampleSA);
	}
}

int main(int argc, char **argv)
{
	checkPopcnt();

	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'a':
				opt::alphabet = arg.str();
				arg.clear(ios::eofbit);
				break;
			case 'c': opt::toStdout = true; break;
			case 'd': opt::decompress = true; break;
			case 's': arg >> opt::sampleSA; break;
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

		ifstream in(fmPath.c_str());
		assert_good(in, fmPath);
		FMIndex fmIndex;
		in >> fmIndex;
		assert_good(in, fmPath);
		in.close();

		ofstream fout;
		if (!opt::toStdout)
			fout.open(faPath.c_str());
		ostream& out = opt::toStdout ? cout : fout;
		assert_good(out, faPath);
		fmIndex.decompress(
				ostream_iterator<FMIndex::value_type>(out, ""));
		out.flush();
		assert_good(out, faPath);

		in.open((faPath + ".fai").c_str());
		FastaIndex faIndex;
		if (in) {
			in >> faIndex;
			faIndex.writeFASTAHeaders(out);
		}
		return 0;
	}


	string path = argv[optind];
	FMIndex fm;
	if (opt::bwt2fm) {
		buildFMIndexFromBWT(fm, path);
	} else {
		if (opt::indexes & opt::FAI)
			indexFasta(path, path + ".fai");

		if ((opt::indexes & opt::FM) == 0)
			return 0;

		buildFMIndex(fm, path);
	}

	if (opt::verbose > 0) {
		size_t n = fm.size();
		ssize_t bytes = getMemoryUsage();
		cerr << "Read " << toSI(n) << "B. "
			"Used " << toSI(bytes) << "B of memory and "
				<< setprecision(3) << (float)bytes / n << " B/bp.\n";
	}

	string fmPath = opt::toStdout ? "-" : path + ".fm";
	cerr << "Writing `" << fmPath << "'...\n";
	ofstream fout;
	if (!opt::toStdout)
		fout.open(fmPath.c_str());
	ostream& out = opt::toStdout ? cout : fout;
	assert_good(out, fmPath);
	out << fm;
	out.flush();
	assert_good(out, fmPath);

	return 0;
}
