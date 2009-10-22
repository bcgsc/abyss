#include "ContigPath.h"
#include "config.h"
#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <map>

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
"  -o, --out=FILE        write result to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static int verbose;
	static string out;
}

static const char* shortopts = "o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "out",         required_argument, NULL, 'o' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
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

static PathMap loadPaths(istream& pathStream)
{
	assert(pathStream.good());
	string line;
	PathMap pathMap;

	while (getline(pathStream, line)) {
		LinearNumKey pathID;
		LinearNumKey contig;
		char dir;
		ContigPath contigPath;
		PathStruct path;
		istringstream s(line);

		int first = -1, last = -1;
		s >> pathID;
		while (s >> contig >> dir) {
			assert(dir == '-' || dir == '+');
			bool isRC = dir == '-';
			MergeNode node;
			node.id = contig;
			node.isRC = isRC;
			contigPath.appendNode(node);
			if (first < 0)
				first = contig;
			last = contig;
		}

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

	refNode = refPath.getNode(refIndex);

	if (currPathStruct.isKeyFirst)
		overlap.secondIsRC = false;
	else {
		overlap.secondIsRC = true;
		currPath.reverse(true);
	}

	currNode = currPath.getNode(currIndex);
	assert(refNode.id == currNode.id);
	overlap.firstIsRC = refNode.isRC != currNode.isRC;

	if (overlap.firstIsRC) {
		refPath.reverse(true);
		refIndex = refPath.getNumNodes() - refIndex - 1;
	}

	while (currIndex < currPath.getNumNodes() &&
			refIndex < refPath.getNumNodes()) {
		if (refPath.getNode(refIndex) == currPath.getNode(currIndex)) {
			currIndex++;
			refIndex++;
		} else
			return;
	}
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

		for (unsigned i = 0; i < path.path.getNumNodes(); i++) {
			PathMapPair result = pathMap.equal_range(path.path.getNode(i).id);
			unsigned dist = distance(result.first, result.second);

			//the first and last element of the path is guarenteed to
			//have at least one entry in the PathMap.
			dist = i == 0 || i == path.path.getNumNodes() - 1 ?
					dist - 1 : dist;

			if (dist > 0)
				for (PathMap::const_iterator currIt = result.first;
						currIt != result.second; currIt++)
					addOverlap(path, currIt->second, i, overlaps);
		}
	}
	return overlaps;
}

static void removeContigs(TrimPathStruct& pathStruct, bool fromBack,
		unsigned overlap)
{
	if (pathStruct.numRemoved[fromBack] >= overlap)
		return;
	ContigPath newPath;
	int max, min;
	if (fromBack) {
		//subtract the found overlap from the original length
		max = pathStruct.path.getNumNodes() - overlap +
			pathStruct.numRemoved[fromBack];
		min = 0;
	} else {
		max = pathStruct.path.getNumNodes();
		min = overlap - pathStruct.numRemoved[fromBack];
	}
	assert(min >= 0 && max <= (int)pathStruct.path.getNumNodes());

	pathStruct.path = pathStruct.path.extractNodes(min, max);
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

static string toString(const ContigPath& path, char sep)
{
	size_t numNodes = path.getNumNodes();
	assert(numNodes > 0);
	MergeNode root = path.getNode(0);
	ostringstream s;
	s << root.id << (root.isRC ? '-' : '+');
	for (size_t i = 1; i < numNodes; ++i) {
		MergeNode mn = path.getNode(i);
		s << sep << mn.id << (mn.isRC ? '-' : '+');
	}
	return s.str();
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
		key = path.path.getNode(0).id;
		pathMap.insert(pair<LinearNumKey, PathStruct>(key, path));
		path.isKeyFirst = false;
		key = path.path.getNode(path.path.getNumNodes() - 1).id;
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
			case 'o': arg >> opt::out; break;
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
		cerr << "There were  " << overlaps.size() / 2 << " overlaps found.\n";
		trimPaths = makeTrimPathMap(pathMap);
		trimOverlaps(trimPaths, overlaps);
		pathMap = makePathMap(trimPaths);
		overlaps = findOverlaps(pathMap);
		trimIterations++;
	}

	for (TrimPathMap::const_iterator trimPathsIt = trimPaths.begin();
			trimPathsIt != trimPaths.end(); trimPathsIt++) {
		if (trimPathsIt->second.path.getNumNodes() < 2) continue;
		string pathString = toString(trimPathsIt->second.path, ' ');
		cout << pathString << '\n';
	}
	cerr << PROGRAM " completed after " << trimIterations
		<< " trim iterations.\n";
}
