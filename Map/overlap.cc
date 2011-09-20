#include "BitUtils.h"
#include "DataLayer/Options.h"
#include "FMIndex.h"
#include "FastaIndex.h"
#include "FastaReader.h"
#include "IOUtil.h"
#include "MemoryUtil.h"
#include "SAM.h"
#include "StringUtil.h"
#include "Uncompress.h"
#include <cassert>
#include <cctype> // for toupper
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "abyss-overlap"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2011 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... FILE\n"
"Find overlaps of [m,k) bases. Sequences are read from FILE.\n"
"Output is written to standard output. The index files FILE.fai\n"
"and FILE.fm will be used if present.\n"
"\n"
"  -m, --min=N             find matches at least N bp [30]\n"
"  -k, --max=N             find matches less than N bp [inf]\n"
"  -j, --threads=N         use N parallel threads [1]\n"
"  -s, --sample=N          sample the suffix array [1]\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** Find matches at least k bp. */
	static unsigned minOverlap = 30;

	/** Find matches less than k bp. */
	static unsigned maxOverlap = UINT_MAX;

	/** Sample the suffix array. */
	static unsigned sampleSA;

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** Verbose output. */
	static int verbose;
}

static const char shortopts[] = "j:k:m:s:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "help", no_argument, NULL, OPT_HELP },
	{ "max", required_argument, NULL, 'k' },
	{ "min", required_argument, NULL, 'm' },
	{ "sample", required_argument, NULL, 's' },
	{ "threads", required_argument, NULL, 'j' },
	{ "verbose", no_argument, NULL, 'v' },
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
		assert(m.qend == qlength);
		a.cigar = ss.str();
	}
	a.mrnm = "*";
	a.mpos = -1;
	a.isize = 0;
	return a;
}

/** Return the mapping of the specified sequence. */
static void findOverlaps(
		const FastaIndex& faIndex, const FMIndex& fmIndex,
		const string& id, const string& seq, bool rc)
{
	size_t pos = seq.size() > opt::maxOverlap
		? seq.size() - opt::maxOverlap + 1 : 1;
	string suffix(seq, pos);
	typedef vector<FMInterval> Matches;
	vector<FMInterval> matches;
	fmIndex.findOverlap(suffix, back_inserter(matches),
			opt::minOverlap);

	for (Matches::const_iterator it = matches.begin();
			it != matches.end(); ++it) {
		const FMInterval& fmi = *it;
		Match m(fmi.qstart + pos, fmi.qend + pos, 0, fmi.u - fmi.l);
		for (unsigned i = fmi.l; i < fmi.u; ++i) {
			m.tstart = fmIndex.locate(i) + 1;
			SAMRecord sam = toSAM(faIndex, m, rc, seq.size());
			sam.qname = id;
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
	}
}

static void findOverlaps(
		const FastaIndex& faIndex, const FMIndex& fmIndex,
		const FastqRecord& rec)
{
	findOverlaps(faIndex, fmIndex, rec.id, rec.seq, false);
	string rcseq(reverseComplement(rec.seq));
	findOverlaps(faIndex, fmIndex, rec.id, rcseq, true);
}

/** Map the sequences of the specified file. */
static void findOverlaps(
		const FastaIndex& faIndex, const FMIndex& fmIndex,
		FastaReader& in)
{
#pragma omp parallel
	for (FastqRecord rec;;) {
		bool good;
#pragma omp critical(in)
		good = in >> rec;
		if (good)
			findOverlaps(faIndex, fmIndex, rec);
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

/** Return the size of the specified file. */
static streampos fileSize(const string& path)
{
	std::ifstream in(path.c_str());
	assert_good(in, path);
	in.seekg(0, std::ios::end);
	assert_good(in, path);
	return in.tellg();
}

/** Check that the indexes are up to date. */
static void checkIndexes(const string& path,
		const FMIndex& fmIndex, const FastaIndex& faIndex)
{
	size_t fastaFileSize = fileSize(path);
	if (fmIndex.size() != fastaFileSize) {
		cerr << PROGRAM ": `" << path << "': "
			"The size of the FM-index, "
			<< fmIndex.size()
			<< " B, does not match the size of the FASTA file, "
			<< fastaFileSize << " B. The index is likely stale.\n";
		exit(EXIT_FAILURE);
	}
	if (faIndex.fileSize() != fastaFileSize) {
		cerr << PROGRAM ": `" << path << "': "
			"The size of the FASTA index, "
			<< faIndex.fileSize()
			<< " B, does not match the size of the FASTA file, "
			<< fastaFileSize << " B. The index is likely stale.\n";
		exit(EXIT_FAILURE);
	}
}

int main(int argc, char** argv)
{
	checkPopcnt();

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
		  case '?':
			die = true; break;
		  case 'j':
			arg >> opt::threads; assert(arg.eof()); break;
		  case 'm':
			arg >> opt::minOverlap; assert(arg.eof()); break;
		  case 'k':
			arg >> opt::maxOverlap; assert(arg.eof()); break;
		  case 's':
			arg >> opt::sampleSA; assert(arg.eof()); break;
		  case 'v':
			opt::verbose++; break;
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
	} else if (argc - optind > 1) {
		cerr << PROGRAM ": too many arguments\n";
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

	assert(opt::minOverlap < opt::maxOverlap);

	const char* fastaFile(argv[--argc]);
	ostringstream ss;
	ss << fastaFile << ".fm";
	string fmPath(ss.str());
	ss.str("");
	ss << fastaFile << ".fai";
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
			cerr << "Reading `" << fastaFile << "'...\n";
		faIndex.index(fastaFile);
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
		buildFMIndex(fmIndex, fastaFile);
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

	// Check that the indexes are up to date.
	checkIndexes(fastaFile, fmIndex, faIndex);

	opt::chastityFilter = false;
	opt::trimMasked = false;
	FastaReader fa(fastaFile, FastaReader::FOLD_CASE);
	findOverlaps(faIndex, fmIndex, fa);

	if (opt::verbose > 0) {
		size_t unique = g_count.unique;
		size_t mapped = unique + g_count.multimapped;
		size_t total = mapped + g_count.unmapped;
		cerr << "Mapped " << mapped << " of " << total << " queries ("
			<< (float)100 * mapped / total << "%)\n"
			<< "Mapped " << unique << " of " << total
			<< " queries uniquely (" << (float)100 * unique / total
			<< "%)\n";
	}

	cout.flush();
	assert_good(cout, "stdout");
	return 0;
}
