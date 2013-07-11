#include "Histogram.h"
#include "IOUtil.h"
#include "MemoryUtil.h"
#include "SAM.h"
#include "StringUtil.h"
#include "Uncompress.h"
#include "UnorderedMap.h"
#include "ContigID.h"
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

using namespace std;

#define PROGRAM "abyss-fixmate"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2013 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [FILE]...\n"
"Write read pairs that map to the same contig to the file SAME.\n"
"Write read pairs that map to different contigs to stdout.\n"
"Alignments may be in FILE(s) or standard input.\n"
"\n"
" Options:\n"
"\n"
"      --no-qname        set the qname to * [default]\n"
"      --qname           do not alter the qname\n"
"  -l, --min-align=N     the minimal alignment size [1]\n"
"  -s, --same=SAME       write properly-paired reads to this file\n"
"  -h, --hist=FILE       write the fragment size histogram to FILE\n"
"  -c, --cov=FILE        write the physical coverage to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static string fragPath;
	static string histPath;
	static string covPath;
	static int qname;
	static int verbose;
}

static const char shortopts[] = "h:c:l:s:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "qname",     no_argument,       &opt::qname, 1 },
	{ "no-qname",  no_argument,       &opt::qname, 0 },
	{ "min-align", required_argument, NULL, 'l' },
	{ "hist",      required_argument, NULL, 'h' },
	{ "cov",       required_argument, NULL, 'c' },
	{ "same",      required_argument, NULL, 's' },
	{ "verbose",   no_argument,       NULL, 'v' },
	{ "help",      no_argument,       NULL, OPT_HELP },
	{ "version",   no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

static struct {
	size_t alignments;
	size_t bothUnaligned;
	size_t oneUnaligned;
	size_t numDifferent;
	size_t numFF;
} stats;

static ofstream g_fragFile;
static Histogram g_histogram;
static ofstream g_covFile;
static vector< vector<int> > g_contigCov;

static void handlePair(SAMRecord& a0, SAMRecord& a1)
{
	if ((a0.isRead1() && a1.isRead1())
			|| (a0.isRead2() && a1.isRead2())) {
		cerr << "error: duplicate read ID `" << a0.qname
			<< (a0.isRead1() ? "/1" : "")
			<< (a0.isRead2() ? "/2" : "")
			<< "'\n";
		exit(EXIT_FAILURE);
	}

	if (!opt::qname)
		a0.qname = a1.qname = "*";

	fixMate(a0, a1);
	if (a0.isUnmapped() && a1.isUnmapped()) {
		// Both reads are unaligned.
		stats.bothUnaligned++;
	} else if (a0.isUnmapped() || a1.isUnmapped()) {
		// One read is unaligned.
		stats.oneUnaligned++;
	} else if (a0.rname != a1.rname) {
		// Different targets.
		stats.numDifferent++;
		// Set the mapping quality of both reads to their minimum.
		a0.mapq = a1.mapq = min(a0.mapq, a1.mapq);
		if (!opt::histPath.empty())
			cout << a0 << '\n' << a1 << '\n';
	} else if (a0.isReverse() == a1.isReverse()) {
		// Same target, FF orientation.
		stats.numFF++;
	} else {
		// Same target, FR or RF orientation.
		SAMRecord a = a0.isize >= 0 ? a0 : a1;
		unsigned inx = get(g_contigNames, a.rname);
		g_contigCov[inx][a.pos]++;
		assert(a.isize >= 0);
		unsigned end = a.pos + a.isize;
		if (end >= g_contigCov[inx].size())
			end = g_contigCov[inx].size() - 1;
		g_contigCov[inx][end]--;
		g_histogram.insert(a0.isReverse() ? a1.isize : a0.isize);
		if (!opt::fragPath.empty()) {
			g_fragFile << a0 << '\n' << a1 << '\n';
			assert(g_fragFile.good());
		}
	}

	if (opt::histPath.empty()) {
		cout << a0 << '\n' << a1 << '\n';
		assert(cout.good());
	}
}

#if SAM_SEQ_QUAL
typedef unordered_map<string, SAMRecord> Alignments;
#else
typedef unordered_map<string, SAMAlignment> Alignments;
#endif

static void printProgress(const Alignments& map)
{
	if (opt::verbose == 0)
		return;

	static size_t prevBuckets;
	if (prevBuckets == 0)
		prevBuckets = map.bucket_count();

	size_t buckets = map.bucket_count();
	if (stats.alignments % 1000000 == 0 || buckets != prevBuckets) {
		prevBuckets = buckets;
		size_t size = map.size();
		cerr << "Read " << stats.alignments << " alignments. "
			<< "Hash load: " << size << " / " << buckets
			<< " = " << (float)size / buckets
			<< " using " << toSI(getMemoryUsage()) << "B." << endl;
	}
}

static void handleAlignment(SAMRecord& sam, Alignments& map)
{
	pair<Alignments::iterator, bool> it = map.insert(
			make_pair(sam.qname, sam));
	if (!it.second) {
#if SAM_SEQ_QUAL
		SAMRecord& a0 = it.first->second;
#else
		SAMRecord a0(it.first->second, it.first->first);
#endif
		handlePair(a0, sam);
		map.erase(it.first);
	}
	stats.alignments++;
	printProgress(map);
}

static void assert_eof(istream& in)
{
	if (in.eof())
		return;
	in.clear();
	string line;
	getline(in, line);
	cerr << "error: `" << line << "'\n";
	exit(EXIT_FAILURE);
}

/** Print physical coverage in wiggle format. */
static void printCov(string file)
{
	ofstream out(file.c_str());
	for (unsigned i = 0; i < g_contigCov.size(); i++) {
		out << "variableStep\tchrom=" << get(g_contigNames, i) << '\n';
		int prev = g_contigCov[i][0];
		if (prev != 0) out << "0\t" << prev << '\n';
		for (unsigned j = 1; j < g_contigCov[i].size(); j++) {
			prev += g_contigCov[i][j];
			if (prev != 0) out << j << "\t" << prev << '\n';
		}
	}
}

void parseTag(string line) {
	stringstream ss(line);
	string tag;
	char type[2];
	string id;
	char delim;
	unsigned length;
	ss >> tag;
	if (tag != "@SQ")
		return;
	while (ss >> type[0] >> type[1] >> delim) {
		assert(delim == ':');
		if (type[0] == 'S' && type[1] == 'N')
			ss >> id;
		else if (type[0] == 'L' && type[1] == 'N')
			ss >> length;
	}

	assert(length > 0);
	assert(id.size() > 0);
	put(g_contigNames, g_contigCov.size(), id);
	g_contigCov.push_back(vector<int>(length));	
}

static void readAlignments(istream& in, Alignments* pMap)
{
	for (SAMRecord sam; in >> ws;) {
		if (in.peek() == '@') {
			string line;
			getline(in, line);
			assert(in);

			parseTag(line);

			cout << line << '\n';
			if (!opt::fragPath.empty())
				g_fragFile << line << '\n';
		} else if (in >> sam)
			handleAlignment(sam, *pMap);
	}
	printCov(opt::covPath);
	assert_eof(in);
}

static void readAlignmentsFile(string path, Alignments* pMap)
{
	if (opt::verbose > 0)
		cerr << "Reading `" << path << "'..." << endl;
	ifstream fin(path.c_str());
	assert_good(fin, path);
	readAlignments(fin, pMap);
	fin.close();
}

/** Return the specified number formatted as a percent. */
static string percent(size_t x, size_t n)
{
	ostringstream ss;
	ss << setw((int)log10(n) + 1) << x;
	if (x > 0)
		ss << "  " << setprecision(3) << (float)100*x/n << '%';
	return ss.str();
}

/** Print statistics of the specified histogram. */
static void printHistogramStats(const Histogram& h)
{
	cerr << "Stats mean: " << setprecision(4) << h.mean() << " "
		"median: " << setprecision(4) << h.median() << " "
		"sd: " << setprecision(4) << h.sd() << " "
		"n: " << h.size() << " "
		"min: " << h.minimum() << " "
		"max: " << h.maximum() << '\n'
		<< h.barplot() << endl;
}

int main(int argc, char* const* argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'l':
				arg >> opt::minAlign;
				break;
			case 's': arg >> opt::fragPath; break;
			case 'h': arg >> opt::histPath; break;
			case 'c': arg >> opt::covPath; break;
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

	if (!opt::fragPath.empty()) {
		g_fragFile.open(opt::fragPath.c_str());
		assert(g_fragFile.is_open());
	}

	Alignments alignments(1);
	if (optind < argc) {
		for_each(argv + optind, argv + argc,
				bind2nd(ptr_fun(readAlignmentsFile), &alignments));
	} else {
		if (opt::verbose > 0)
			cerr << "Reading from standard input..." << endl;
		readAlignments(cin, &alignments);
	}
	if (opt::verbose > 0)
		cerr << "Read " << stats.alignments << " alignments" << endl;

	// Print the unpaired alignments.
	if (opt::histPath.empty()) {
		for (Alignments::iterator it = alignments.begin();
				it != alignments.end(); it++) {
#if SAM_SEQ_QUAL
			SAMRecord& a0 = it->second;
#else
			SAMRecord a0(it->second, it->first);
#endif
			cout << a0 << '\n';
			assert(cout.good());
		}
	}

	unsigned numRF = g_histogram.count(INT_MIN, 0);
	unsigned numFR = g_histogram.count(1, INT_MAX);
	size_t sum = alignments.size()
		+ stats.bothUnaligned + stats.oneUnaligned
		+ numFR + numRF + stats.numFF
		+ stats.numDifferent;
	cerr <<
		"Mateless   " << percent(alignments.size(), sum) << "\n"
		"Unaligned  " << percent(stats.bothUnaligned, sum) << "\n"
		"Singleton  " << percent(stats.oneUnaligned, sum) << "\n"
		"FR         " << percent(numFR, sum) << "\n"
		"RF         " << percent(numRF, sum) << "\n"
		"FF         " << percent(stats.numFF, sum) << "\n"
		"Different  " << percent(stats.numDifferent, sum) << "\n"
		"Total      " << sum << endl;
	
	if (alignments.size() == sum) {
		cerr << PROGRAM ": error: All reads are mateless. This "
			"can happen when first and second read IDs do not match."
			<< endl;
		exit(EXIT_FAILURE);
	}

	if (!opt::fragPath.empty())
		g_fragFile.close();

	if (!opt::histPath.empty()) {
		ofstream histFile(opt::histPath.c_str());
		assert(histFile.is_open());
		histFile << g_histogram;
		assert(histFile.good());
		histFile.close();
	}

	if (opt::verbose > 0) {
		size_t numTotal = numFR + numRF;

		// Print the statistics of the forward-reverse distribution.
		if ((float)numFR / numTotal > 0.001) {
			Histogram h = g_histogram;
			h.eraseNegative();
			h.removeNoise();
			h.removeOutliers();
			cerr << "FR ";
			printHistogramStats(h.trimFraction(0.0001));
		}

		// Print the statistics of the reverse-forward distribution.
		if ((float)numRF / numTotal > 0.001) {
			Histogram h = g_histogram.negate();
			h.eraseNegative();
			h.removeNoise();
			h.removeOutliers();
			cerr << "RF ";
			printHistogramStats(h.trimFraction(0.0001));
		}
	}

	if (stats.numFF > numFR && stats.numFF > numRF) {
		cerr << "error: The mate pairs of this library are oriented "
			"forward-forward (FF), which is not supported by ABySS."
			<< endl;
		exit(EXIT_FAILURE);
	}

	return 0;
}
