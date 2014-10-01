#include "BitUtil.h"
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
#include <boost/algorithm/string/join.hpp>
#include <boost/tuple/tuple.hpp>
#include <algorithm>
#include <cassert>
#include <cctype> // for toupper
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <stdint.h>
#include <utility>
#include <queue>
#if _OPENMP
# include <omp.h>
#endif
#if _SQL
#include "DataBase/Options.h"
#include "DataBase/DB.h"
#endif

using namespace std;
using namespace boost;
using namespace boost::algorithm;

#define PROGRAM "abyss-map"

#if _SQL
DB db;
#endif

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2014 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... QUERY... TARGET\n"
"Map the sequences of the files QUERY to those of the file TARGET.\n"
"The index files TARGET.fai and TARGET.fm will be used if present.\n"
"\n"
" Options:\n"
"\n"
"  -l, --min-align=N       find matches at least N bp [1]\n"
"  -j, --threads=N         use N parallel threads [1]\n"
"  -s, --sample=N          sample the suffix array [1]\n"
"  -d, --dup               identify and print duplicate sequence\n"
"                          IDs between QUERY and TARGET\n"
"      --order             print alignments in the same order as\n"
"                          read from QUERY\n"
"      --no-order          print alignments ASAP [default]\n"
"      --multi             Align unaligned segments of primary\n"
"                          alignment\n"
"      --no-multi          don't Align unaligned segments [default]\n"
"      --SS                expect contigs to be oriented correctly\n"
"      --no-SS             no assumption about contig orientation\n"
"      --rc                map the sequence and its reverse complement [default]\n"
"      --no-rc             do not map the reverse complement sequence\n"
"      --chastity          discard unchaste reads\n"
"      --no-chastity       do not discard unchaste reads [default]\n"
"  -v, --verbose           display verbose output\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
#if _SQL
"      --db=FILE           specify path of database repository in FILE\n"
"      --library=NAME      specify library NAME for database\n"
"      --strain=NAME       specify strain NAME for database\n"
"      --species=NAME      specify species NAME for database\n"
#endif
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
#if _SQL
	string url;
	dbVars metaVars;
#endif
	/** Find matches at least k bp. */
	static unsigned k;

	/** Sample the suffix array. */
	static unsigned sampleSA;

	/** The number of parallel threads. */
	static unsigned threads = 1;

	/** Run a strand-specific RNA-Seq alignments. */
	static int ss;

	/** Do not map the sequence's reverse complement. */
	static int norc;

	/** Identify duplicate and subsumed sequences. */
	static bool dup = false;

	/** Align unaligned segments of primary alignment. */
	static int multi;

	/** Ensure output order matches input order. */
	static int order;

	/** Verbose output. */
	static int verbose;
}

static const char shortopts[] = "j:k:l:s:dv";

#if _SQL
enum { OPT_HELP = 1, OPT_VERSION, OPT_DB, OPT_LIBRARY, OPT_STRAIN, OPT_SPECIES };
#else
enum { OPT_HELP = 1, OPT_VERSION };
#endif

static const struct option longopts[] = {
	{ "sample", required_argument, NULL, 's' },
	{ "min-align", required_argument, NULL, 'l' },
	{ "dup", no_argument, NULL, 'd' },
	{ "threads", required_argument, NULL, 'j' },
	{ "order", no_argument, &opt::order, 1 },
	{ "no-order", no_argument, &opt::order, 0 },
	{ "multi", no_argument, &opt::multi, 1 },
	{ "no-multi", no_argument, &opt::multi, 0 },
	{ "SS", no_argument, &opt::ss, 1 },
	{ "no-SS", no_argument, &opt::ss, 0 },
	{ "rc", no_argument, &opt::norc, 0 },
	{ "no-rc", no_argument, &opt::norc, 1 },
	{ "verbose", no_argument, NULL, 'v' },
	{ "chastity", no_argument, &opt::chastityFilter, 1 },
	{ "no-chastity", no_argument, &opt::chastityFilter, 0 },
	{ "help", no_argument, NULL, OPT_HELP },
	{ "version", no_argument, NULL, OPT_VERSION },
#if _SQL
	{ "db", required_argument, NULL, OPT_DB },
	{ "library", required_argument, NULL, OPT_LIBRARY },
	{ "strain", required_argument, NULL, OPT_STRAIN },
	{ "species", required_argument, NULL, OPT_SPECIES },
#endif
	{ NULL, 0, NULL, 0 }
};

/** Counts. */
static struct {
	unsigned unique;
	unsigned multimapped;
	unsigned unmapped;
	unsigned suboptimal;
	unsigned subunmapped;
} g_count;

typedef FMIndex::Match Match;

#if SAM_SEQ_QUAL
static string toXA(const FastaIndex& faIndex,
		const FMIndex& fmIndex, const Match& m, bool rc,
		unsigned qlength, unsigned seq_start)
{
	if (m.size() == 0)
		return "";
	FastaIndex::SeqPos seqPos = faIndex[fmIndex[m.l]];
	string rname = seqPos.get<0>().id;
	int pos = seqPos.get<1>() + 1;

	// Set the mapq to the alignment score.
	assert(m.qstart < m.qend);
	unsigned matches = m.qend - m.qstart;

	unsigned qstart = seq_start + m.qstart;
	unsigned qend = m.qend + seq_start;
	unsigned short flag = rc ? SAMAlignment::FREVERSE : 0;

	ostringstream ss;
	if (qstart > 0)
		ss << qstart << 'S';
	ss << matches << 'M';
	if (qend < qlength)
		ss << qlength - qend << 'S';
	string cigar = ss.str();

	stringstream xa_tag;
	xa_tag << rname << ',' << pos << ',' << cigar << ",0," << flag;

	return xa_tag.str();
}
#endif

/** Return a SAM record of the specified match. */
static SAMRecord toSAM(const FastaIndex& faIndex,
		const FMIndex& fmIndex, const Match& m, bool rc,
		unsigned qlength)
{
	SAMRecord a;
	if (m.size() == 0) {
		// No hit.
		a.rname = "*";
		a.pos = -1;
		a.flag = SAMAlignment::FUNMAP;
		a.mapq = 0;
		a.cigar = "*";
	} else {
		FastaIndex::SeqPos seqPos = faIndex[fmIndex[m.l]];
		a.rname = seqPos.get<0>().id;
		a.pos = seqPos.get<1>();
		a.flag = rc ? SAMAlignment::FREVERSE : 0;

		// Set the mapq to the alignment score.
		assert(m.qstart < m.qend);
		unsigned matches = m.qend - m.qstart;
		assert (m.num != 0);
		a.mapq = m.size() > 1 || m.num > 1 ? 0 : min(matches, 254U);

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

/** Return the position of the current contig. */
static size_t getMyPos(const Match& m, const FastaIndex& faIndex,
		const FMIndex& fmIndex, const string& id)
{
	for (size_t i = m.l; i < m.u; i++) {
		if (faIndex[fmIndex[i]].get<0>().id == id)
			return fmIndex[i];
	}
	return fmIndex[m.l];
}

/** Return the earlies position of all contigs in m. */
static size_t getMinPos(const Match& m, size_t maxLen,
		const FastaIndex& faIndex, const FMIndex& fmIndex)
{
	size_t minPos = numeric_limits<size_t>::max();
	for (size_t i = m.l; i < m.u; i++) {
		size_t pos = fmIndex[i];
		if (faIndex[pos].get<0>().size == maxLen && pos < minPos)
			minPos = fmIndex[i];
	}
	return minPos;
}

/** Return the largest length of all contig in m. */
static size_t getMaxLen(const Match& m, const FastaIndex& faIndex,
		const FMIndex& fmIndex)
{
	size_t maxLen = 0;
	for (size_t i = m.l; i < m.u; i++) {
		size_t len = faIndex[fmIndex[i]].get<0>().size;
		if (len > maxLen)
			maxLen = len;
	}
	return maxLen;
}

/** Print the current contig id if it is not the lartest and earliest
 * contig in m. */
static void printDuplicates(const Match& m, const Match& rcm,
		const FastaIndex& faIndex, const FMIndex& fmIndex,
		const FastqRecord& rec)
{
	size_t myLen = m.qspan();
	size_t maxLen;
	if (opt::ss || opt::norc)
		maxLen = getMaxLen(m, faIndex, fmIndex);
	else
		maxLen = max(getMaxLen(m, faIndex, fmIndex),
				getMaxLen(rcm, faIndex, fmIndex));
	if (myLen < maxLen) {
#pragma omp atomic
		g_count.multimapped++;
#pragma omp critical(cout)
		{
			cout << rec.id << '\n';
			assert_good(cout, "stdout");
		}
		return;
	}
	size_t myPos = getMyPos(m, faIndex, fmIndex, rec.id);
	size_t minPos;
	if (opt::ss || opt::norc)
		minPos = getMinPos(m, maxLen, faIndex, fmIndex);
	else
		minPos = min(getMinPos(m, maxLen, faIndex, fmIndex),
				getMinPos(rcm, maxLen, faIndex, fmIndex));
	if (myPos > minPos) {
#pragma omp atomic
		g_count.multimapped++;
#pragma omp critical(cout)
		{
			cout << rec.id << '\n';
			assert_good(cout, "stdout");
		}
	}
#pragma omp atomic
	g_count.unique++;
	return;
}

pair<Match, Match> findMatch(const FMIndex& fmIndex, const string& seq)
{
	Match m = fmIndex.find(seq,
			opt::dup ? seq.length() : opt::k);
	if (opt::norc)
		return make_pair(m, Match());

	string rcqseq = reverseComplement(seq);
	Match rcm;
	if (opt::ss)
		rcm = fmIndex.find(rcqseq,
				opt::dup ? rcqseq.length() : opt::k);
	else
		rcm = fmIndex.find(rcqseq,
				opt::dup ? rcqseq.length() : m.qspan());
	return make_pair(m, rcm);
}

static queue<string> g_pq;

/** Return the mapping of the specified sequence. */
static void find(const FastaIndex& faIndex, const FMIndex& fmIndex,
		const FastqRecord& rec)
{
	if (rec.seq.empty()) {
		cerr << PROGRAM ": error: "
			"the sequence `" << rec.id << "' is empty\n";
		exit(EXIT_FAILURE);
	}

	Match m, rcm;
	tie(m, rcm) = findMatch(fmIndex, rec.seq);

	if (opt::dup) {
		printDuplicates(m, rcm, faIndex, fmIndex, rec);
		return;
	}

	bool rc;
	if (opt::ss) {
		rc = rec.id.size() > 2
			&& rec.id.substr(rec.id.size()-2) == "/1";
		bool prc = rcm.qspan() > m.qspan();
		if (prc != rc && ((rc && rcm.size() > 0)
					|| (!rc && m.size() > 0)))
#pragma omp atomic
			g_count.suboptimal++;
		if (prc != rc && ((rc && rcm.size() == 0 && m.size() > 0)
				|| (!rc && m.size() == 0 && rcm.size() > 0)))
#pragma omp atomic
			g_count.subunmapped++;
	} else {
		rc = rcm.qspan() > m.qspan();

		// if both matches are the same length, sum up the number of times
		// each were seen.
		if (rcm.qspan() == m.qspan())
			rc ? rcm.num += m.num : m.num += rcm.num;
	}

	vector<string> alts;
	Match mm = rc ? rcm : m;
	string mseq = rc ? reverseComplement(rec.seq) : rec.seq;

#if SAM_SEQ_QUAL
	if (opt::multi) {
		if (mm.qstart > 0) {
			string seq = mseq.substr(0, mm.qstart);
			Match m1, rcm1;
			tie(m1, rcm1) = findMatch(fmIndex, seq);
			bool rc1 = rcm1.qspan() > m1.qspan();
			string xa = toXA(faIndex, fmIndex, rc1 ? rcm1 : m1,
					rc ^ rc1, mseq.size(), 0);
			if (xa != "")
				alts.push_back(xa);
		}

		if (mm.qend < mseq.size()) {
			string seq =
				mseq.substr(mm.qend, mseq.length() - mm.qend);
			Match m2, rcm2;
			tie(m2, rcm2) = findMatch(fmIndex, seq);
			bool rc2 = rcm2.qspan() > m2.qspan();
			string xa = toXA(faIndex, fmIndex, rc2 ? rcm2 : m2,
					rc ^ rc2, mseq.size(), mm.qend);
			if (xa != "")
				alts.push_back(xa);
		}
	}
#endif

	SAMRecord sam = toSAM(faIndex, fmIndex, mm, rc,
			rec.seq.size());
	if (rec.id[0] == '@') {
		cerr << PROGRAM ": error: "
			"the query ID `" << rec.id << "' is invalid since it "
			"begins with `@'\n";
		exit(EXIT_FAILURE);
	}
	sam.qname = rec.id;

#if SAM_SEQ_QUAL
	sam.seq = mseq;
	sam.qual = rec.qual.empty() ? "*" : rec.qual;
	if (rc)
		reverse(sam.qual.begin(), sam.qual.end());
#endif

	bool print = opt::order == 0;
	do {
#pragma omp critical(cout)
		{
#pragma omp critical(g_pq)
			if (!print) {
				print = g_pq.front() == rec.id;
				if (print)
					g_pq.pop();
			}
			if (print) {
				cout << sam;
#if SAM_SEQ_QUAL
				if (alts.size() > 0)
					cout << "\tXA:Z:" << join(alts, ";");
#endif
				cout << '\n';
				assert_good(cout, "stdout");
			}
		}
	} while (!print); // spinlock :(

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
		{
			good = in >> rec;
			if (opt::order) {
#pragma omp critical(g_pq)
				g_pq.push(rec.id);
			}
		}
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
	string commandLine;
	{
		ostringstream ss;
		char** last = argv + argc - 1;
		copy(argv, last, ostream_iterator<const char *>(ss, " "));
		ss << *last;
		commandLine = ss.str();
	}

	opt::chastityFilter = false;
	opt::trimMasked = false;

#if _SQL
	opt::metaVars.resize(3);
#endif

	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'j': arg >> opt::threads; break;
			case 'k': case 'l':
				arg >> opt::k;
				break;
			case 's': arg >> opt::sampleSA; break;
			case 'd': opt::dup = true; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
#if _SQL
			case OPT_DB:
				arg >> opt::url; break;
			case OPT_LIBRARY:
				arg >> opt::metaVars[0]; break;
			case OPT_STRAIN:
				arg >> opt::metaVars[1]; break;
			case OPT_SPECIES:
				arg >> opt::metaVars[2]; break;
#endif
		}
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

#ifndef SAM_SEQ_QUAL
# define SAM_SEQ_QUAL 0
#endif

	if (opt::multi && !SAM_SEQ_QUAL) {
		cerr << PROGRAM ": multiple alignments not supported with "
			"this install. Recompile ABySS with `./configure "
			"--enable-samseqqual'.\n";
		die = true;
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

#if _SQL
	init(db,
			opt::url,
			opt::verbose,
			PROGRAM,
			opt::getCommand(argc, argv),
			opt::metaVars
	);
	addToDb(db, "K", opt::k);
	addToDb(db, "SS", opt::ss);
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
	if (opt::verbose > 0) {
		ssize_t bytes = getMemoryUsage();
		if (bytes > 0)
			cerr << "Using " << toSI(bytes) << "B of memory and "
				<< setprecision(3) << (float)bytes / faIndex.size()
				<< " B/sequence.\n";
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
#if _SQL
	addToDb(db, "readContigs", faIndex.size());
#endif

	// Check that the indexes are up to date.
	checkIndexes(targetFile, fmIndex, faIndex);

	if (!opt::dup) {
		// Write the SAM header.
		cout << "@HD\tVN:1.4\n"
			"@PG\tID:" PROGRAM "\tPN:" PROGRAM "\tVN:" VERSION "\t"
			"CL:" << commandLine << '\n';
		faIndex.writeSAMHeader(cout);
		cout.flush();
		assert_good(cout, "stdout");
	} else if (opt::verbose > 0)
		cerr << "Identifying duplicates.\n";

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
#if _SQL
		addToDb(db, "read_alignments_initial", total);
		addToDb(db, "mapped", mapped);
		addToDb(db, "mapped_uniq", unique);
#endif
		if (opt::ss) {
			cerr << "Mapped " << g_count.suboptimal
				<< " (" << (float)100 * g_count.suboptimal / total << "%)"
				<< " reads to the opposite strand of the optimal mapping.\n"
				<< "Made " << g_count.subunmapped << " ("
				<< (float)100 * g_count.subunmapped / total << "%)"
				<< " unmapped suboptimal decisions.\n";
#if _SQL
			addToDb(db, "reads_map_ss", g_count.suboptimal);
#endif
		} else {
#if _SQL
			addToDb(db, "reads_map_ss", 0);
#endif
		}
	}

	cout.flush();
	assert_good(cout, "stdout");
	return 0;
}
