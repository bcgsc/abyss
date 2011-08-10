#include "DataLayer/Options.h"
#include "FMIndex.h"
#include "FastaIndex.h"
#include "FastaInterleave.h"
#include "FastaReader.h"
#include "IOUtil.h"
#include "SAM.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string>
#include <utility>
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "abyss-map"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... QUERY... TARGET\n"
"Map the sequences of the files QUERY to those of the file TARGET.\n"
"The TARGET file must be indexed. The files TARGET.fai and\n"
"TARGET.fm are required.\n"
"\n"
"  -k, --score=N           find matches at least N bp [1]\n"
"  -j, --threads=N       use N parallel threads [1]\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** Find matches at least k bp. */
	static unsigned k;

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** Verbose output. */
	static int verbose;
}

static const char shortopts[] = "j:k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "score", required_argument, NULL, 'k' },
	{ "threads", required_argument, NULL, 'j' },
	{ "verbose", no_argument, NULL, 'v' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Counts. */
static struct {
	unsigned unique;
	unsigned multimapped;
	unsigned unmapped;
} g_count;

/** Return a SAM record of the specified match. */
static SAMRecord toSAM(const FastaIndex& faIndex,
		const Match& m, bool rc, unsigned qlength)
{
	SAMRecord a;
	if (m.count == 0) {
		// No hit.
		a.rname = "*";
		a.pos = -1;
		a.flag = SAMAlignment::FUNMAP;
		a.mapq = 0;
		a.cigar = "*";
	} else {
		pair<string, size_t> idPos = faIndex[m.tstart];
		a.rname = idPos.first;
		a.pos = idPos.second;
		a.flag = rc ? SAMAlignment::FREVERSE : 0;
		a.mapq = m.count > 1 ? 0 : 255;
		ostringstream ss;
		if (m.qstart > 0)
			ss << m.qstart << 'S';
		ss << m.qend - m.qstart << 'M';
		if (m.qend < qlength)
			ss << qlength - m.qend << 'S';
		a.cigar = ss.str();
	}
	a.mrnm = "*";
	a.mpos = -1;
	a.isize = 0;
	return a;
}

/** Return the mapping of the specified sequence. */
static void find(const FastaIndex& faIndex, const FMIndex& fmIndex,
		const FastqRecord& rec)
{
	Match m = fmIndex.find(rec.seq, opt::k);
	string rcqseq = reverseComplement(rec.seq);
	Match rcm = fmIndex.find(rcqseq, opt::k);
	bool rc = rcm.qspan() > m.qspan();

	SAMRecord sam = toSAM(faIndex, rc ? rcm : m, rc, rec.seq.size());
	sam.qname = rec.id;
#if SAM_SEQ_QUAL
	sam.seq = rc ? rcqseq : rec.seq;
	sam.qual = rec.qual.empty() ? "*" : rec.qual;
	if (rc)
		reverse(sam.qual.begin(), sam.qual.end());
#endif
#pragma omp critical(cout)
	cout << sam << '\n';

	if (sam.isUnmapped())
#pragma omp atomic
		g_count.unmapped++;
	else if (sam.mapq == 0)
#pragma omp atomic
		g_count.multimapped++;
	else
#pragma omp atomic
		g_count.unique++;
}

/** Map the sequences of the specified file. */
static void find(const FastaIndex& faIndex, const FMIndex& fmIndex,
		FastaInterleave& in)
{
#pragma omp parallel
	for (FastqRecord rec;;) {
		bool good;
#pragma omp critical(in)
		good = in >> rec;
		if (good)
			find(faIndex, fmIndex, rec);
		else
			break;
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
			case 'j': arg >> opt::threads; assert(arg.eof()); break;
			case 'k': arg >> opt::k; assert(arg.eof()); break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

#if _OPENMP
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

	const char* targetFile(argv[--argc]);
	ostringstream ss;
	ss << targetFile << ".fm";
	string fmPath(ss.str());
	ss.str("");
	ss << targetFile << ".fai";
	string faiPath(ss.str());

	ifstream in;

	// Read the FASTA index.
	FastaIndex faIndex;
	in.open(faiPath.c_str());
	if (in) {
		if (opt::verbose > 0)
			cerr << "Reading `" << faiPath << "'...\n";
		in >> faIndex;
		assert(in.eof());
		in.close();
	} else {
		if (opt::verbose > 0)
			cerr << "Reading `" << targetFile << "'...\n";
		faIndex.index(targetFile);
	}

	// Read the FM index.
	FMIndex fmIndex;
	in.open(fmPath.c_str());
	if (in) {
		if (opt::verbose > 0)
			cerr << "Reading `" << fmPath << "'...\n";
		assert_good(in, fmPath);
		in >> fmIndex;
		assert_good(in, fmPath);
		in.close();
	} else {
		fmIndex.setAlphabet("\nACGT");
		fmIndex.buildIndex(targetFile);
	}

	opt::chastityFilter = false;
	opt::trimMasked = false;
	FastaInterleave fa(argv + optind, argv + argc,
			FastaReader::FOLD_CASE);
	find(faIndex, fmIndex, fa);

	if (opt::verbose > 0) {
		size_t unique = g_count.unique;
		size_t mapped = unique + g_count.multimapped;
		size_t total = mapped + g_count.unmapped;
		cerr << "Mapped " << mapped << " of " << total << " reads ("
			<< (float)100 * mapped / total << "%)\n"
			<< "Mapped " << unique << " of " << total
			<< " reads uniquely (" << (float)100 * unique / total
			<< "%)\n";
	}

	cout.flush();
	assert_good(cout, "stdout");
	return 0;
}
