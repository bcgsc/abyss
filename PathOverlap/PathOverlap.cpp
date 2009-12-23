#include "config.h"
#include "ContigPath.h"
#include "PairUtils.h"
#include <algorithm>
#include <cassert>
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
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... PATH\n"
"\n"
"  -r, --repeats=FILE write repeat contigs to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static int verbose;
	static string repeatContigs;
}

static const char* shortopts = "r:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "repeats",      required_argument, NULL, 'r' },
	{ "verbose",      no_argument,       NULL, 'v' },
	{ "help",         no_argument,       NULL, OPT_HELP },
	{ "version",      no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

struct PathStruct {
	LinearNumKey pathID;
	ContigPath path;
	bool isKeyFirst;
};

struct OverlapStruct {
	LinearNumKey firstID, secondID;
	bool firstIsRC, secondIsRC;
	unsigned overlap;
};

struct TrimPathStruct {
	ContigPath path;
	vector<unsigned> numRemoved;
};

typedef multimap<LinearNumKey, PathStruct> PathMap;
typedef pair<PathMap::iterator, PathMap::iterator> PathMapPair;
typedef map<LinearNumKey, TrimPathStruct> TrimPathMap;
typedef vector<OverlapStruct> OverlapVec;

static set<LinearNumKey> s_trimmedContigs;

static PathMap loadPaths(istream& pathStream)
{
	assert(pathStream.good());
	string line;
	PathMap pathMap;

	while (getline(pathStream, line)) {
		istringstream s(line);
		LinearNumKey pathID;
		s >> pathID;

		ContigPath contigPath;
		copy(istream_iterator<MergeNode>(s),
				istream_iterator<MergeNode>(),
				back_inserter(contigPath));
		assert(s.eof());
		LinearNumKey first = contigPath.front().id();
		LinearNumKey last = contigPath.back().id();

		PathStruct path;
		path.path = contigPath;
		path.pathID = pathID;
		path.isKeyFirst = true;
		pathMap.insert(pair<LinearNumKey, PathStruct>(first, path));
		path.isKeyFirst = false;
		pathMap.insert(pair<LinearNumKey, PathStruct>(last, path));
	}
	return pathMap;
}

static void addOverlap(const PathStruct& refPathStruct,
		const PathStruct& currPathStruct, unsigned refIndex,
		OverlapVec& overlaps)
{
	if (refPathStruct.pathID == currPathStruct.pathID) return;

	OverlapStruct overlap;
	ContigPath refPath = refPathStruct.path;
	ContigPath currPath = currPathStruct.path;
	MergeNode refNode, currNode;
	unsigned currIndex = 0;

	refNode = refPath[refIndex];

	if (currPathStruct.isKeyFirst)
		overlap.secondIsRC = false;
	else {
		overlap.secondIsRC = true;
		currPath.reverseComplement();
	}

	currNode = currPath[currIndex];
	assert(refNode.id() == currNode.id());
	overlap.firstIsRC = refNode.sense() != currNode.sense();

	if (overlap.firstIsRC) {
		refPath.reverseComplement();
		refIndex = refPath.size() - refIndex - 1;
	}

//	bool isRepeat = currIndex == 0 && refIndex == 0;

	while (currIndex < currPath.size() &&
			refIndex < refPath.size()) {
		if (refPath[refIndex] == currPath[currIndex]) {
			currIndex++;
			refIndex++;
		} else {
//			if (isRepeat)
//				break;
//			else
				return;
		}
	}
	//overlap.firstID = overlap.firstID != isRepeat;
	overlap.overlap = currIndex;
	overlap.firstID = refPathStruct.pathID;
	overlap.secondID = currPathStruct.pathID;
	overlaps.push_back(overlap);
}

static OverlapVec findOverlaps(PathMap& pathMap)
{
	OverlapVec overlaps;
	for (PathMap::const_iterator pathIt = pathMap.begin();
			pathIt != pathMap.end(); pathIt++) {
		const PathStruct& path = pathIt->second;
		if (!path.isKeyFirst)
			continue;

		for (unsigned i = 0; i < path.path.size(); i++) {
			PathMapPair result = pathMap.equal_range(
					path.path[i].id());
			unsigned dist = distance(result.first, result.second);

			//the first and last element of the path is guarenteed to
			//have at least one entry in the PathMap.
			dist = i == 0 || i == path.path.size() - 1 ?
					dist - 1 : dist;

			if (dist > 0)
				for (PathMap::const_iterator currIt = result.first;
						currIt != result.second; currIt++)
					addOverlap(path, currIt->second, i, overlaps);
		}
	}
	return overlaps;
}

static void addContigsToSet(ContigPath& contigs, int min, int max)
{
	for (int x = 0; x < min; x++)
		s_trimmedContigs.insert(contigs[x].id());

	for (int x = max; (unsigned)x < contigs.size(); x++)
		s_trimmedContigs.insert(contigs[x].id());
}

static void removeContigs(TrimPathStruct& pathStruct, bool fromBack,
		unsigned overlap)
{
	if (pathStruct.numRemoved[fromBack] >= overlap ||
			pathStruct.path.size() == 0)
		return;
	ContigPath newPath;
	int max, min;
	if (fromBack) {
		//subtract the found overlap from the original length
		max = pathStruct.path.size() - overlap +
			pathStruct.numRemoved[fromBack];
		min = 0;
		//assert(max > 0);
	} else {
		max = pathStruct.path.size();
		min = overlap - pathStruct.numRemoved[fromBack];
	}
	assert(min >= 0 && max <= (int)pathStruct.path.size());

	if (min > max) return;
	addContigsToSet(pathStruct.path, min, max);
	pathStruct.path = ContigPath(&pathStruct.path[min],
			&pathStruct.path[max]);
	pathStruct.numRemoved[fromBack] = overlap;
}

static void trimOverlaps(TrimPathMap& pathMap, OverlapVec& overlaps)
{
	for (OverlapVec::const_iterator overlapIt = overlaps.begin();
			overlapIt != overlaps.end(); overlapIt++) {
		TrimPathStruct& firstPath =
			pathMap.find(overlapIt->firstID)->second;
		removeContigs(firstPath, !overlapIt->firstIsRC, overlapIt->overlap);
		TrimPathStruct& secondPath =
			pathMap.find(overlapIt->secondID)->second;
		removeContigs(secondPath, overlapIt->secondIsRC, overlapIt->overlap);
	}
}

static TrimPathMap makeTrimPathMap(const PathMap& pathMap)
{
	TrimPathMap trimMap;
	for (PathMap::const_iterator pathIt = pathMap.begin();
			pathIt != pathMap.end(); pathIt++) {
		if (pathIt->second.isKeyFirst) {
			TrimPathStruct path;
			path.path = pathIt->second.path;
			vector<unsigned> numRemoved(2);
			path.numRemoved = numRemoved;
			trimMap[pathIt->second.pathID] = path;
		}
	}
	return trimMap;
}

static PathMap makePathMap(const TrimPathMap& trimPathMap)
{
	PathMap pathMap;
	for (TrimPathMap::const_iterator trimIt = trimPathMap.begin();
			trimIt != trimPathMap.end(); trimIt++) {
		PathStruct path;
		LinearNumKey key;
		path.path = trimIt->second.path;
		path.pathID = trimIt->first;
		path.isKeyFirst = true;
		if (path.path.size() < 1)
			continue;
		key = path.path[0].id();
		pathMap.insert(pair<LinearNumKey, PathStruct>(key, path));
		path.isKeyFirst = false;
		key = path.path[path.path.size() - 1].id();
		pathMap.insert(pair<LinearNumKey, PathStruct>(key, path));
	}
	return pathMap;
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

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	PathMap pathMap;
	if (optind < argc) {
		ifstream fin(argv[argc - 1]);
		pathMap = loadPaths(fin);
	} else {
		pathMap = loadPaths(cin);
	}

	int trimIterations = 0;
	OverlapVec overlaps;
	TrimPathMap trimPaths;

	overlaps = findOverlaps(pathMap);
	while (overlaps.size() > 0) {
		cerr << "There were " << overlaps.size() / 2 << " overlaps found.\n";
		trimPaths = makeTrimPathMap(pathMap);
		trimOverlaps(trimPaths, overlaps);
		pathMap = makePathMap(trimPaths);
		overlaps = findOverlaps(pathMap);
		trimIterations++;
	}

	for (TrimPathMap::const_iterator trimPathsIt = trimPaths.begin();
			trimPathsIt != trimPaths.end(); trimPathsIt++) {
		if (trimPathsIt->second.path.size() < 2) continue;
		cout << trimPathsIt->first << '\t'
			<< trimPathsIt->second.path << '\n';
	}

	if (!opt::repeatContigs.empty()) {
		ofstream out(opt::repeatContigs.c_str());
		for (set<LinearNumKey>::const_iterator it
				= s_trimmedContigs.begin();
				it != s_trimmedContigs.end(); ++it)
			out << g_contigIDs.key(*it) << '\n';
		assert(out.good());
	}

	cerr << PROGRAM " completed after " << trimIterations
		<< " trim iterations.\n";
}
