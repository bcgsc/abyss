#include "config.h"
#include "ContigID.h"
#include "ContigLength.h"
#include "ContigPath.h"
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstring> // for strerror
#include <cstdlib>
#include <numeric>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <getopt.h>
#include <map>
#include <set>

using namespace std;

#define PROGRAM "PathOverlap"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Tony Raymond.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... LEN PATH\n"
"Find paths that overlap\n"
"  LEN   lengths of the contigs\n"
"  PATH  sequences of contig IDs\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -r, --repeats=FILE    write repeat contigs to FILE\n"
"      --dot             output overlaps in dot format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by readContigLengths
	
	/** Output overlaps in dot format. Do not perform trimming. */
	static int dot;

	/** Output the IDs of contigs in overlaps to this file. */
	static string repeatContigs;

	static int verbose;
}

static const char* shortopts = "k:r:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",         required_argument, NULL, 'k' },
	{ "dot",          no_argument,       &opt::dot, 1, },
	{ "repeats",      required_argument, NULL, 'r' },
	{ "verbose",      no_argument,       NULL, 'v' },
	{ "help",         no_argument,       NULL, OPT_HELP },
	{ "version",      no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** Lengths of contigs in k-mer. */
static vector<unsigned> g_contigLengths;

/** Return the length of this contig in k-mer. */
unsigned ContigNode::length() const
{
	return ambiguous() ? m_id : g_contigLengths.at(id());
}

/** An alignment seed. */
struct Seed {
	LinearNumKey pathID;
	const ContigPath& path;
	bool isKeyFirst;
	Seed(LinearNumKey id, const ContigPath& path, bool front)
		: pathID(id), path(path), isKeyFirst(front) { }
};

/** An alignment result. */
struct Overlap {
	LinearNumKey firstID, secondID;
	bool firstIsRC, secondIsRC;

	/** Overlap measured in number of contigs. */
	unsigned short overlap;

	/** Overlap measured in bp. */
	int distance;

	friend ostream& operator <<(ostream& out, Overlap o)
	{
		return out <<
			'"' << o.firstID << (o.firstIsRC ? '-' : '+') << "\" -> "
			"\"" << o.secondID << (o.secondIsRC ? '-' : '+') << "\""
			" [d=" << o.distance << "]";
	}
};

/** A path and its maximum overlap with other paths. */
struct Path {
	ContigPath path;
	unsigned numRemoved[2];
	Path(const ContigPath& path) : path(path)
	{
		numRemoved[0] = numRemoved[1] = 0;
	}
};

/** The contig IDs that have been removed from paths. */
static set<LinearNumKey> s_trimmedContigs;

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

typedef map<LinearNumKey, Path> Paths;

/** Read contig paths from the specified stream. */
static Paths readPaths(const string& inPath)
{
	ifstream fin(inPath.c_str());
	if (inPath != "-")
		assert_open(fin, inPath);
	istream& in = inPath == "-" ? cin : fin;

	assert(in.good());
	Paths paths;
	LinearNumKey id;
	ContigPath path;
	while (in >> id >> path) {
		bool inserted = paths.insert(make_pair(id, path)).second;
		assert(inserted);
		(void)inserted;
	}
	assert(in.eof());
	return paths;
}

typedef multimap<ContigNode, Seed> SeedMap;

/** Index the first and last contig of each path to facilitate finding
 * overlaps between paths. */
static SeedMap makeSeedMap(const Paths& paths)
{
	SeedMap seedMap;
	for (Paths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		const ContigPath& path = it->second.path;
		if (path.empty())
			continue;
		seedMap.insert(make_pair(path.front(),
					Seed(it->first, path, true)));
		seedMap.insert(make_pair(~path.back(),
					Seed(it->first, path, false)));
	}
	return seedMap;
}

/** Add the number of k-mer in two contigs. */
static unsigned addLength(unsigned addend, const ContigNode& contig)
{
	return addend + contig.length();
}

/** Check whether path starts with the sequence [first, last). */
static bool startsWith(ContigPath path, bool rc,
		ContigPath::const_iterator first,
		ContigPath::const_iterator last)
{
	if (rc)
		path.reverseComplement();
	assert(*first == path.front());
	assert(first < last);
	return unsigned(last - first) > path.size() ? false
		: equal(first, last, path.begin());
}

/** Check whether path starts with the sequence [first, last). */
static unsigned findOverlap(ContigPath::const_iterator first,
		ContigPath::const_iterator last,
		ContigPath path, bool rc, int &distance)
{
	if (!startsWith(path, rc, first, last))
		return 0;
	distance = -accumulate(first, last, opt::k-1, addLength);
	return last - first;
}

typedef vector<Overlap> Overlaps;

/** Find every path that overlaps with the specified path. */
static void findOverlaps(const SeedMap& seedMap,
		LinearNumKey id, bool rc, const ContigPath& path,
		Overlaps& overlaps)
{
	for (ContigPath::const_iterator it = path.begin();
			it != path.end(); ++it) {
		if (it->ambiguous())
			continue;

		pair<SeedMap::const_iterator, SeedMap::const_iterator>
			range = seedMap.equal_range(*it);
		for (SeedMap::const_iterator seed = range.first;
				seed != range.second; ++seed) {
			if (id == seed->second.pathID)
				continue;
			Overlap o;
			o.secondIsRC = !seed->second.isKeyFirst;
			o.overlap = findOverlap(it, path.end(),
					   seed->second.path, o.secondIsRC, o.distance);
			if (o.overlap > 0) {
				o.firstID = id;
				o.firstIsRC = rc;
				o.secondID = seed->second.pathID;
				overlaps.push_back(o);
			}
		}
	}
}

/** Find every pair of overlapping paths. */
static Overlaps findOverlaps(const Paths& paths)
{
	SeedMap seedMap = makeSeedMap(paths);

	Overlaps overlaps;
	for (Paths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		const ContigPath& path = it->second.path;
		findOverlaps(seedMap, it->first, false, path, overlaps);

		ContigPath rc = path;
		rc.reverseComplement();
		findOverlaps(seedMap, it->first, true, rc, overlaps);
	}
	return overlaps;
}

/** Find the largest overlap for each contig. */
static void determineMaxOverlap(Path& pathStruct,
		bool fromBack, unsigned overlap)
{
	unsigned& curOverlap = pathStruct.numRemoved[fromBack];
	curOverlap = max(curOverlap, overlap);
}

/** Record the trimmed contigs. */
static void recordTrimmedContigs(
		ContigPath::const_iterator first,
		ContigPath::const_iterator last)
{
	for (ContigPath::const_iterator it = first; it != last; ++it)
		if (!it->ambiguous())
			s_trimmedContigs.insert(it->id());
}

/** Remove the overlapping portion of the specified contig. */
static void removeContigs(Path& o)
{
	assert(o.numRemoved[0] <= o.path.size());
	assert(o.numRemoved[1] <= o.path.size());
	unsigned first = o.numRemoved[0];
	unsigned last = o.path.size() - o.numRemoved[1];
	if (first < last) {
		recordTrimmedContigs(o.path.begin(), o.path.begin() + first);
		recordTrimmedContigs(o.path.begin() + last, o.path.end());
		o.path.erase(o.path.begin() + last, o.path.end());
		o.path.erase(o.path.begin(), o.path.begin() + first);
	} else {
		recordTrimmedContigs(o.path.begin(), o.path.end());
		o.path.clear();
	}
}

/** Find the largest overlap for each contig and remove it. */
static void trimOverlaps(Paths& paths, const Overlaps& overlaps)
{
	for (Paths::iterator it = paths.begin(); it != paths.end(); ++it)
		it->second.numRemoved[0] = it->second.numRemoved[1] = 0;

	for (Overlaps::const_iterator it = overlaps.begin();
			it != overlaps.end(); ++it) {
		determineMaxOverlap(paths.find(it->firstID)->second,
				!it->firstIsRC, it->overlap);
		determineMaxOverlap(paths.find(it->secondID)->second,
				it->secondIsRC, it->overlap);
	}

	for (Paths::iterator it = paths.begin(); it != paths.end(); ++it)
		removeContigs(it->second);
}

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
			case 'r': arg >> opt::repeatContigs; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				cout << USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				cout << VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
		}
	}

	if (opt::k <= 0) {
		cerr << PROGRAM ": missing -k,--kmer option\n";
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

	g_contigLengths = readContigLengths(argv[optind++]);
	Paths paths = readPaths(argv[optind++]);

	if (opt::dot) {
		Overlaps overlaps = findOverlaps(paths);
		cout << "digraph path_overlap {\n";
		copy(overlaps.begin(), overlaps.end(),
				ostream_iterator<Overlap>(cout, "\n"));
		cout << "}\n";
		return 0;
	}

	for (Overlaps overlaps = findOverlaps(paths);
			!overlaps.empty(); overlaps = findOverlaps(paths)) {
		cerr << "Found " << overlaps.size() / 2 << " overlaps.\n";
		trimOverlaps(paths, overlaps);
	}

	for (Paths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		if (it->second.path.size() < 2)
			continue;
		cout << it->first << '\t' << it->second.path << '\n';
	}

	if (!opt::repeatContigs.empty()) {
		ofstream out(opt::repeatContigs.c_str());
		for (set<LinearNumKey>::const_iterator it
				= s_trimmedContigs.begin();
				it != s_trimmedContigs.end(); ++it)
			out << ContigID(*it) << '\n';
		assert(out.good());
	}

	return 0;
}
