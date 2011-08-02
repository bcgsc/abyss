#include "DataLayer/Options.h"
#include "FMIndex.h"
#include "FastaIndex.h"
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
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** Find matches at least k bp. */
	static unsigned k;

	static int verbose;
}

static const char shortopts[] = "k:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "score", required_argument, NULL, 'k' },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

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
	assert(m.qstart <= m.qend);
	assert(m.qend <= rec.seq.size());

	bool rc = false;
	string rcqseq;
	if (m.count == 0) {
		rcqseq = reverseComplement(rec.seq);
		m = fmIndex.find(rcqseq, opt::k);
		rc = m.count > 0;
	}

	SAMRecord sam = toSAM(faIndex, m, rc, rec.seq.size());
	sam.qname = rec.id;
#if SAM_SEQ_QUAL
	sam.seq = rc ? rcqseq : rec.seq;
	sam.qual = rec.qual.empty() ? "*" : rec.qual;
	if (rc)
		reverse(sam.qual.begin(), sam.qual.end());
#endif
	cout << sam << '\n';
}

/** Map the sequences of the specified file. */
static void find(const FastaIndex& faIndex, const FMIndex& fmIndex,
		const char* path)
{
	FastaReader in(path, FastaReader::FOLD_CASE);
	for (FastqRecord rec; in >> rec;)
		find(faIndex, fmIndex, rec);
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

	const char* targetFile(argv[--argc]);
	ostringstream ss;
	ss << targetFile << ".fm";
	FMIndex fmIndex(ss.str());

	ss.str("");
	ss << targetFile << ".fai";
	FastaIndex faIndex(ss.str());

	opt::chastityFilter = false;
	opt::trimMasked = false;
	for (char** p = argv + optind; p != argv + argc; ++p)
		find(faIndex, fmIndex, *p);

	cout.flush();
	assert_good(cout, "stdout");
	return 0;
}
