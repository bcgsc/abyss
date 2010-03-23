#include "config.h"
#include "ContigPath.h"
#include "PairUtils.h"
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstring> // for strerror
#include <cstdlib>
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
"Usage: " PROGRAM " [OPTION]... PATH\n"
"\n"
"  -r, --repeats=FILE    write repeat contigs to FILE\n"
"      --dot             output overlaps in dot format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	/** Output overlaps in dot format. Do not perform trimming. */
	static int dot;

	/** Output the IDs of contigs in overlaps to this file. */
	static string repeatContigs;

	static int verbose;
}

static const char* shortopts = "r:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "dot",          no_argument,       &opt::dot, 1, },
	{ "repeats",      required_argument, NULL, 'r' },
	{ "verbose",      no_argument,       NULL, 'v' },
	{ "help",         no_argument,       NULL, OPT_HELP },
	{ "version",      no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

struct PathStruct {
	LinearNumKey pathID;
	const ContigPath& path;
	bool isKeyFirst;
	PathStruct(LinearNumKey id, const ContigPath& path, bool front)
		: pathID(id), path(path), isKeyFirst(front) { }
};

struct OverlapStruct {
	LinearNumKey firstID, secondID;
	bool firstIsRC, secondIsRC;
	unsigned overlap;

	friend ostream& operator <<(ostream& out, OverlapStruct o)
	{
		return out <<
			'"' << o.firstID << (o.firstIsRC ? '-' : '+') << "\" -> "
			"\"" << o.secondID << (o.secondIsRC ? '-' : '+') << "\""
			" [label = " << o.overlap << "];";
	}
};

struct TrimPathStruct {
	ContigPath path;
	unsigned numRemoved[2];
	TrimPathStruct(const ContigPath& path) : path(path)
	{
		numRemoved[0] = numRemoved[1] = 0;
	}
};

typedef map<LinearNumKey, TrimPathStruct> TrimPathMap;
typedef vector<OverlapStruct> OverlapVec;

/** The contig IDs that have been removed from paths. */
static set<LinearNumKey> s_trimmedContigs;

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

/** Read contig paths from the specified stream. */
static TrimPathMap readPaths(const string& inPath)
{
	ifstream fin(inPath.c_str());
	if (inPath != "-")
		assert_open(fin, inPath);
	istream& in = inPath == "-" ? cin : fin;

	assert(in.good());
	TrimPathMap paths;
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

typedef multimap<ContigNode, PathStruct> PathMap;

/** Index the first and last contig of each path to facilitate finding
 * overlaps between paths. */
static PathMap makePathMap(const TrimPathMap& paths)
{
	PathMap pathMap;
	for (TrimPathMap::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		const ContigPath& path = it->second.path;
		if (path.empty())
			continue;
		pathMap.insert(make_pair(path.front(),
					PathStruct(it->first, path, true)));
		pathMap.insert(make_pair(~path.back(),
					PathStruct(it->first, path, false)));
	}
	return pathMap;
}

/** Check whether these two paths overlaps. */
static unsigned findOverlap(ContigPath::const_iterator first,
		ContigPath::const_iterator last,
		ContigPath path2, bool rc)
{
	if (rc)
		path2.reverseComplement();
	assert(*first == path2.front());
	unsigned size1 = last - first;
	bool overlap = size1 < path2.size()
		? equal(first, last, path2.begin())
		: equal(path2.begin(), path2.end(), first);
	return overlap ? min(size1, path2.size()) : 0;
}

/** Find every path that overlaps with the specified path. */
static void findOverlaps(const PathMap& pathMap,
		LinearNumKey id, bool rc, const ContigPath& path,
		OverlapVec& overlaps)
{
	for (ContigPath::const_iterator it = path.begin();
			it != path.end(); ++it) {
		if (it->ambiguous())
			continue;

		pair<PathMap::const_iterator, PathMap::const_iterator>
			range = pathMap.equal_range(*it);
		for (PathMap::const_iterator mapIt = range.first;
				mapIt != range.second; ++mapIt) {
			if (id == mapIt->second.pathID)
				continue;
			OverlapStruct o;
			o.secondIsRC = !mapIt->second.isKeyFirst;
			o.overlap = findOverlap(it, path.end(),
					   mapIt->second.path, o.secondIsRC);
			if (o.overlap > 0) {
				o.firstID = id;
				o.firstIsRC = rc;
				o.secondID = mapIt->second.pathID;
				overlaps.push_back(o);
			}
		}
	}
}

/** Find every pair of overlapping paths. */
static OverlapVec findOverlaps(const TrimPathMap& paths)
{
	PathMap pathMap = makePathMap(paths);

	OverlapVec overlaps;
	for (TrimPathMap::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		const ContigPath& path = it->second.path;
		findOverlaps(pathMap, it->first, false, path, overlaps);

		ContigPath rc = path;
		rc.reverseComplement();
		findOverlaps(pathMap, it->first, true, rc, overlaps);
	}
	return overlaps;
}

/** Find the largest overlap for each contig. */
static void determineMaxOverlap(TrimPathStruct& pathStruct,
		bool fromBack, unsigned overlap)
{
	unsigned& curOverlap = pathStruct.numRemoved[fromBack];
	curOverlap = max(curOverlap, overlap);
}

/** Remove the overlapping portion of the specified contig. */
static void removeContigs(TrimPathStruct& o)
{
	assert(o.numRemoved[0] <= o.path.size());
	assert(o.numRemoved[1] <= o.path.size());
	unsigned first = o.numRemoved[0];
	unsigned last = o.path.size() - o.numRemoved[1];
	if (first < last) {
		transform(o.path.begin(), o.path.begin() + first,
				inserter(s_trimmedContigs, s_trimmedContigs.begin()),
				mem_fun_ref(&ContigNode::id));
		transform(o.path.begin() + last, o.path.end(),
				inserter(s_trimmedContigs, s_trimmedContigs.begin()),
				mem_fun_ref(&ContigNode::id));
		o.path.erase(o.path.begin() + last, o.path.end());
		o.path.erase(o.path.begin(), o.path.begin() + first);
	} else {
		transform(o.path.begin(), o.path.end(),
				inserter(s_trimmedContigs, s_trimmedContigs.begin()),
				mem_fun_ref(&ContigNode::id));
		o.path.clear();
	}
}

/** Find the largest overlap for each contig and remove it. */
static void trimOverlaps(TrimPathMap& paths,
		const OverlapVec& overlaps)
{
	for (TrimPathMap::iterator it = paths.begin();
			it != paths.end(); ++it)
		it->second.numRemoved[0] = it->second.numRemoved[1] = 0;

	for (OverlapVec::const_iterator overlapIt = overlaps.begin();
			overlapIt != overlaps.end(); ++overlapIt) {
		determineMaxOverlap(paths.find(overlapIt->firstID)->second,
				!overlapIt->firstIsRC, overlapIt->overlap);
		determineMaxOverlap(paths.find(overlapIt->secondID)->second,
				overlapIt->secondIsRC, overlapIt->overlap);
	}

	for (TrimPathMap::iterator it = paths.begin();
			it != paths.end(); ++it)
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

	if (argc - optind < 1) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	TrimPathMap paths = readPaths(argv[optind++]);

	if (opt::dot) {
		OverlapVec overlaps = findOverlaps(paths);
		cout << "digraph path_overlap {\n";
		copy(overlaps.begin(), overlaps.end(),
				ostream_iterator<OverlapStruct>(cout, "\n"));
		cout << "}\n";
		return 0;
	}

	for (OverlapVec overlaps = findOverlaps(paths);
			!overlaps.empty(); overlaps = findOverlaps(paths)) {
		cerr << "Found " << overlaps.size() / 2 << " overlaps.\n";
		trimOverlaps(paths, overlaps);
	}

	for (TrimPathMap::const_iterator it = paths.begin();
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
			out << g_contigIDs.key(*it) << '\n';
		assert(out.good());
	}
	
	return 0;
}
