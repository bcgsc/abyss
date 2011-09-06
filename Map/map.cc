#include "DataLayer/Options.h"
#include "FMIndex.h"
#include "FastaIndex.h"
#include "FastaInterleave.h"
#include "FastaReader.h"
#include "IOUtil.h"
#include "MemoryUtil.h"
#include "SAM.h"
#include "StringUtil.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <cctype> // for toupper
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
"The index files TARGET.fai and TARGET.fm will be used if present.\n"
"\n"
"  -k, --score=N           find matches at least N bp [1]\n"
"  -j, --threads=N         use N parallel threads [1]\n"
"  -s, --sample=N          sample the suffix array [1]\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** Find matches at least k bp. */
	static unsigned k;

	/** Sample the suffix array. */
	static unsigned sampleSA;

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** Verbose output. */
	static int verbose;
}

static const char shortopts[] = "j:k:s:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "sample", required_argument, NULL, 's' },
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

		// Set the mapq to the alignment score.
		assert(m.qstart < m.qend);
		unsigned matches = m.qend - m.qstart;
		a.mapq = m.count > 1 ? 0 : min(matches, 255U);

		ostringstream ss;
		if (m.qstart > 0)
			ss << m.qstart << 'S';
		ss << matches << 'M';
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
	Match rcm = fmIndex.find(rcqseq, max(opt::k, m.qspan()));
	bool rc = rcm.qspan() > m.qspan();

	SAMRecord sam = toSAM(faIndex, rc ? rcm : m, rc, rec.seq.size());
	sam.qname = rec.id;
#if SAM_SEQ_QUAL
	sam.seq = rc ? rcqseq : rec.seq;
	sam.qual = rec.qual.empty() ? "*" : rec.qual;
	if (rc)
		reverse(sam.qual.begin(), sam.qual.end());
#endif

	if (m.qstart == rec.seq.size() - rcm.qend
			&& m.qspan() == rcm.qspan()) {
		// This matching sequence maps to both strands.
		sam.mapq = 0;
	}

#pragma omp critical(cout)
	{
		cout << sam << '\n';
		assert_good(cout, "stdout");
	}

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

/** Build an FM index of the specified file. */
static void buildFMIndex(FMIndex& fm, const char* path)
{
	if (opt::verbose > 0)
		std::cerr << "Reading `" << path << "'...\n";
	std::vector<FMIndex::value_type> s;
	readFile(path, s);

	size_t MAX_SIZE = numeric_limits<FMIndex::sais_size_type>::max();
	if (s.size() > MAX_SIZE) {
		std::cerr << PROGRAM << ": `" << path << "', "
			<< toSI(s.size())
			<< "B, must be smaller than " 
			<< toSI(MAX_SIZE) << "B\n";
		exit(EXIT_FAILURE);
	}

	transform(s.begin(), s.end(), s.begin(), ::toupper);
	fm.setAlphabet("-ACGT");
	fm.assign(s.begin(), s.end());
}

int main(int argc, char** argv)
{
	string commandLine;
	{
		ostringstream ss;
		char** last = argv + argc - 1;
		copy(argv, last, ostream_iterator<const char *>(ss, " "));
		ss << *last;
		commandLine = ss.str();
	}

	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'j': arg >> opt::threads; assert(arg.eof()); break;
			case 'k': arg >> opt::k; assert(arg.eof()); break;
			case 's': arg >> opt::sampleSA; assert(arg.eof()); break;
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
	} else
		buildFMIndex(fmIndex, targetFile);
	if (opt::sampleSA > 1)
		fmIndex.sampleSA(opt::sampleSA);

	if (opt::verbose > 0) {
		size_t bp = fmIndex.size();
		cerr << "Read " << toSI(bp) << "B in "
			<< faIndex.size() << " contigs.\n";
		ssize_t bytes = getMemoryUsage();
		if (bytes > 0)
			cerr << "Using " << toSI(bytes) << "B of memory and "
				<< setprecision(3) << (float)bytes / bp << " B/bp.\n";
	}

	// Write the SAM header.
	cout << "@HD\tVN:1.4\n"
		"@PG\tID:" PROGRAM "\tPN:" PROGRAM "\tVN:" VERSION "\t"
		"CL:" << commandLine << '\n';
	faIndex.writeSAMHeader(cout);
	cout.flush();
	assert_good(cout, "stdout");

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
