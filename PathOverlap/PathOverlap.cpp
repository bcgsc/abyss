#include "ContigPath.h"
#include "config.h"
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
};
typedef multimap<LinearNumKey, PathStruct> PathMap;
typedef pair<PathMap::iterator, PathMap::iterator> PathMapPair;

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
		pathMap.insert(pair<LinearNumKey, PathStruct>(first, path));
		pathMap.insert(pair<LinearNumKey, PathStruct>(last, path));
	}
	return pathMap;
}

static void printOverlap(const PathStruct& refPathStruct,
		const PathStruct& currPathStruct, unsigned refIndex)
{
	if (refPathStruct.pathID == currPathStruct.pathID) return;
	ContigPath refPath = refPathStruct.path;
	ContigPath currPath = currPathStruct.path;
	MergeNode refNode, currNode;
	bool refIsRC, currIsRC;
	unsigned currIndex = 0;

	refNode = refPath.getNode(refIndex);
	assert(refPath.getNode(refIndex).id ==
			currPath.getNode(currIndex).id || 
			refPath.getNode(refIndex).id ==
			currPath.getNode(currPath.getNumNodes()-1).id );
	//this may be too simple...
	if (refNode.id == currPath.getNode(0).id)
		currIsRC = false;
	else {
		currIsRC = true;
		currPath.reverse(true);
	}

	currNode = currPath.getNode(currIndex);
	refIsRC = refNode.isRC != currNode.isRC;

	if (refIsRC) {
		refPath.reverse(true);
		refIndex = refPath.getNumNodes() - refIndex - 1;
	}

	bool isMatch = true;
	while (currIndex < currPath.getNumNodes() &&
			refIndex < refPath.getNumNodes()) {
		isMatch = refPath.getNode(refIndex) ==
			currPath.getNode(currIndex);
		if (!isMatch)
			return;
		currIndex++;
		refIndex++;
	}
	cout << refPathStruct.pathID << "," << refIsRC << " "
			<< currPathStruct.pathID << "," << currIsRC << " "
			<< currIndex << '\n';
}

static void findOverlaps(PathMap& pathMap)
{
	for (PathMap::const_iterator pathIt = pathMap.begin();
			pathIt != pathMap.end(); pathIt++) {
		//cout << pathIt->second << '\n';

		const PathStruct& path = pathIt->second;

		for (unsigned i = 0; i < path.path.getNumNodes(); i++) {
			PathMapPair result = pathMap.equal_range(path.path.getNode(i).id);
			unsigned dist = distance(result.first, result.second);

			//the first and last element of the path is guarenteed to
			//have at least one entry in the PathMap.
			dist = i == 0 || i == path.path.getNumNodes() - 1 ?
					dist - 1 : dist;

			if (dist > 0) {
				for (PathMap::const_iterator currIt = result.first;
						currIt != result.second; currIt++) {
					printOverlap(path, currIt->second, i);
				}
			}
		}
	}
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
	if (opt::out.empty())
		pathMap = loadPaths(cin);
	else {
		ifstream fin(opt::out.c_str());
		pathMap = loadPaths(fin);
	}

	findOverlaps(pathMap);
}
