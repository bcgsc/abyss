#include "config.h"
#include "ContigPath.h"
#include "FastaWriter.h"
#include "PackedSeq.h"
#include "PairedAlgorithms.h"
#include "PairUtils.h"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <list>
#include <set>
#include <sstream>
#include <string>

using namespace std;

#define PROGRAM "MergePaths"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2009 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... CONTIG PATH\n"
"Merge paths and contigs.\n"
"  CONTIG  contigs in FASTA format\n"
"  PATH    paths through these contigs\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -o, --out=FILE        write result to FILE\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static int verbose;
	static string out;
}

static const char* shortopts = "k:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "out",         required_argument, NULL, 'o' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

typedef std::list<MergeNode> MergeNodeList;
typedef std::map<LinearNumKey, ContigPath*> ContigPathMap;

// Functions
void readIDIntPair(std::string str, LinearNumKey& id, int& i);
void readPathsFromFile(std::string pathFile, ContigPathMap& contigPathMap);
void parsePathLine(std::string pathLine, LinearNumKey& id, extDirection& dir, ContigPath& path);
void linkPaths(LinearNumKey id, ContigPathMap& contigPathMap);
void mergePath(LinearNumKey cID, const ContigVec& sourceContigs,
		const ContigPath& mergeRecord, int count, int kmer,
		FastaWriter* writer);
void mergeSequences(Sequence& rootContig, const Sequence& otherContig, extDirection dir, bool isReversed, size_t kmer);
bool extractMinCoordSet(LinearNumKey anchor, ContigPath& path, vector<size_t>& coords);
bool checkPathConsistency(LinearNumKey path1Root, LinearNumKey path2Root, ContigPath& path1, ContigPath& path2, size_t& startP1, size_t& endP1, size_t& startP2, size_t& endP2);
void addPathNodesToList(MergeNodeList& list, ContigPath& path);

static bool gDebugPrint;

static set<size_t> getContigIDs(const ContigPathMap& contigPathMap)
{
	set<size_t> seen;
	for (ContigPathMap::const_iterator it = contigPathMap.begin();
			it != contigPathMap.end(); it++) {
		const ContigPath &cp = *it->second;
		size_t nodes = cp.getNumNodes();
		for (size_t i = 0; i < nodes; i++)
			seen.insert(cp.getNode(i).id);
	}
	return seen;
}

int main(int argc, char** argv)
{
	bool die = false;
	for (char c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'k': arg >> opt::k; break;
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

	if (opt::k <= 0) {
		cerr << PROGRAM ": missing -k,--kmer option\n";
		die = true;
	}

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	} else if (argc - optind > 2) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	gDebugPrint = opt::verbose > 1;

	string contigFile(argv[optind++]);
	string pathFile(argv[optind++]);

	// Set up the ID->sequence mapping
	ContigVec contigVec;
	PairedAlgorithms::readContigVec(contigFile, contigVec);
	
	// Read the paths file
	ContigPathMap contigPathMap;
	readPathsFromFile(pathFile, contigPathMap);

	// link the paths together
	for (ContigPathMap::const_iterator iter = contigPathMap.begin();
			iter != contigPathMap.end(); ++iter) {
		linkPaths(iter->first, contigPathMap);

		if (gDebugPrint)
			std::cout << "Pseudo final path from " << iter->first
				<< ' ' << iter->second << " is " << *iter->second << std::endl;
	}

	set<ContigPath*> unique;
	for (ContigPathMap::const_iterator it = contigPathMap.begin();
			it != contigPathMap.end(); ++it)
		unique.insert(it->second);

	FastaWriter writer(opt::out.c_str());

	set<size_t> seen = getContigIDs(contigPathMap);
	for (size_t i = 0; i < contigVec.size(); i++) {
		if (seen.count(i) == 0)
			writer.WriteSequence(contigVec[i].seq,
					i, contigVec[i].coverage);
	}

	int id = contigVec.size();
	for (set<ContigPath*>::const_iterator it = unique.begin();
			it != unique.end(); ++it)
		mergePath((**it).getNode(0).id, contigVec, **it, id++,
				opt::k, &writer);

	return 0;
}

static void assert_open(std::ifstream& f, const std::string& p)
{
	if (f.is_open())
		return;
	std::cerr << p << ": " << strerror(errno) << std::endl;
	exit(EXIT_FAILURE);
}

void readPathsFromFile(std::string pathFile, ContigPathMap& contigPathMap)
{
	std::ifstream pathStream(pathFile.c_str());
	assert_open(pathStream, pathFile);

	while(!pathStream.eof() && pathStream.peek() != EOF)
	{
		// read a line
		std::string pathRecord;
		getline(pathStream, pathRecord);
		
		// parse the line
		LinearNumKey id;
		extDirection dir;
		ContigPath path;
		parsePathLine(pathRecord, id, dir, path);
		//std::cout << "Parsed " << id << " dir " << dir << std::endl;

		if (contigPathMap.find(id) == contigPathMap.end()) {
			MergeNode rootNode = {id, 0};
			(contigPathMap[id] = new ContigPath)->appendNode(rootNode);
		}
		if (dir == 0)
			contigPathMap[id]->appendPath(path);
		else {
			path.reverse(false);
			contigPathMap[id]->prependPath(path);
		}
		//cout << "Adding: " << path << " to: " << id << "," << dir << " to get: " << *contigPathMap[id] << endl;
	}

	pathStream.close();
}

void linkPaths(LinearNumKey id, ContigPathMap& contigPathMap)
{
	ContigPath* refCanonical = contigPathMap[id];

	if(gDebugPrint) std::cout << "Initial canonical path (" << id << ") " << *refCanonical << "\n";

	// Build the initial list of nodes to attempt to merge in
	MergeNodeList mergeInList;
	addPathNodesToList(mergeInList, *refCanonical);

	MergeNodeList::iterator iter = mergeInList.begin();
	while(!mergeInList.empty()) {
		if(iter->id != id) {
			if(gDebugPrint) std::cout << "CHECKING NODE " << iter->id << "(" << iter->isRC << ")\n";

			// Check if the current node to merge has any paths to/from it
			ContigPathMap::iterator findIter = contigPathMap.find(iter->id);
			if (findIter != contigPathMap.end() && refCanonical != findIter->second) {
				// Make the full path of the child node
				ContigPath* childCanonPath = findIter->second;

				if(gDebugPrint) std::cout << " ref: " << refCanonical
					<< ' ' << *refCanonical << "\n";
				if(gDebugPrint) std::cout << "  in: " <<
					childCanonPath << ' ' << *childCanonPath << "\n";

				size_t s1, s2, e1, e2;
				bool validMerge = checkPathConsistency(id, iter->id,
					*refCanonical, *childCanonPath, s1, e1, s2, e2);

				if(validMerge) {
					// Extract the extra nodes from the child path that can be added in
					ContigPath prependNodes = childCanonPath->extractNodes(0, s2);
					ContigPath appendNodes = childCanonPath->extractNodes(e2+1, childCanonPath->getNumNodes());

					// Add the nodes to the list of contigs to try to merge in
					addPathNodesToList(mergeInList, prependNodes);
					addPathNodesToList(mergeInList, appendNodes);

					//std::cout << "PPN " << prependNodes << "\n";
					//std::cout << "APN " << appendNodes << "\n";

					// Add the nodes to the ref contig
					refCanonical->prependPath(prependNodes);
					refCanonical->appendPath(appendNodes);

					if(gDebugPrint) std::cout << " new: " <<
						refCanonical << ' ' << *refCanonical << "\n";

					// Alias all pointers to the child to the reference
					MergeNodeList::iterator mergeIt = mergeInList.begin();
					while (mergeIt != mergeInList.end()) {
						ContigPathMap::iterator delIt = contigPathMap.find(mergeIt->id);
						if (delIt->second == childCanonPath) {
							delIt->second = refCanonical;
							if (mergeIt != iter) {
								mergeInList.erase(mergeIt++);
							} else
								++mergeIt;
						} else
							++mergeIt;
					}
				}
			}
		}

		// Erase the iterator and move forward
		mergeInList.erase(iter++);
	}
}

//
// Check if the two paths are consistent
// They are consistent if there is an identical subpath thats belongs to both nodes and that subpath is terminal wrt to its super path
//
bool checkPathConsistency(LinearNumKey path1Root, LinearNumKey path2Root, ContigPath& path1, ContigPath& path2, size_t& startP1, size_t& endP1, size_t& startP2, size_t& endP2)
{
	(void)path1Root;
	// Find the provisional minimal index set by choosing the closest index pair of the root nodes from each path
	// Since each path must contain each root node, if the range of these indices are different
	// the paths must be different

	assert(path1.getNumNodes() != 0 && path2.getNumNodes() != 1);

	// Extract the minimal coordinates of the root nodes in the paths
	// These coordinates should have the same size
	vector<size_t> coords1, coords2;
	bool valid1 = extractMinCoordSet(path2Root, path1, coords1);
	bool valid2 = extractMinCoordSet(path2Root, path2, coords2);

	// Check that the nodes are both found and the range is the same size
	if(!valid1 || !valid2) //trivially inconsistent
		return false;

	//printf("Init  coords: [%zu-%zu] [%zu-%zu]\n", startP1, endP1, startP2, endP2);
	bool lowValid = true;
	bool highValid = true;
	bool flipped = false;
	size_t max1 = path1.getNumNodes() - 1;
	size_t max2 = path2.getNumNodes() - 1;
	for (unsigned i = 0; i < coords1.size(); i++) {
		for (unsigned j = 0; j < coords2.size(); j++) {
			if(path1.getNode(coords1[i]).isRC != path2.getNode(coords2[j]).isRC) {
				// Flip the path if node direction is different
				path2.reverse(true);
				flipped = !flipped;
			}
			startP1 = coords1[i];
			endP1 = coords1[i];
			if (flipped) {
				startP2 = max2 - coords2[j];
				endP2 = max2 - coords2[j];
			} else {
				startP2 = coords2[j];
				endP2 = coords2[j];
			}

			lowValid = true;
			while(1) {
				if(path1.getNode(startP1).id != path2.getNode(startP2).id) {
					// The nodes no longer match, this path is not valid
					lowValid = false;
					break;
				}

				// Can we expand any further?
				if(startP1 == 0 || startP2 == 0)
					break;

				startP1--;
				startP2--;
			}

			// high coordinates
			highValid = true;
			while(1) {
				if(path1.getNode(endP1).id != path2.getNode(endP2).id) {
					// The nodes no longer match, this path is not valid
					highValid = false;
					break;
				}

				// Can we expand any further?
				if(endP1 == max1 || endP2 == max2)
					break;

				endP1++;
				endP2++;
			}
			if (lowValid && highValid)
				break;
		}
		if (lowValid && highValid)
			break;
	}

	// Check if there was an actual mismatch in the nodes
	if(!lowValid || !highValid) {
		if(gDebugPrint) printf("Invalid path match!\n");
		if(gDebugPrint) std::cout << "Path1 (" << path1Root << ") " << path1 << std::endl;
		if(gDebugPrint) std::cout << "Path2 (" << path2Root << ") " << path2 << std::endl;
		return false;
	}
	
	//printf("Final coords: [%zu-%zu] [%zu-%zu]\n", startP1, endP1, startP2, endP2);
	
	// Sanity assert, at this point one of the low coordniates should be zero and one of the high coordinates should be (size -1)	
	assert(startP1 == 0 || startP2 == 0);
	assert(endP1 == max1 || endP2 == max2);
	
	// Ensure the consistency is met
	size_t count = endP1 - startP1;
	assert(endP2 - startP2 == count);	
	for(size_t c = 0; c < count; ++c)
	{
		if(path1.getNode(startP1 + c).id != path2.getNode(startP2 + c).id)
		{
			if(gDebugPrint) printf("Internal path mismatch\n");
			return false;
		}
	}
	
	// If we got to this point there is a legal subpath that describes both nodes and they can be merged
	return true;
}

// Extract the minimal coordinate set of the indices of (c1, c2) from path.
// Returns true if a valid coordinate set is found, false otherwise
bool extractMinCoordSet(LinearNumKey anchor, ContigPath& path,
		vector<size_t>& coords)
{
	size_t maxIdx = path.getNumNodes();
	for(size_t idx = 0; idx < maxIdx; ++idx) {
		size_t tIdx = maxIdx - idx - 1;
		if(path.getNode(tIdx).id == anchor)
			coords.push_back(tIdx);
	}

	if(coords.empty()) // anchor coord not found
		return false;

	return true;

	/*
	printf("	found %zu %zu %zu %zu\n", coords1[0], coords1[1], coords2[0], coords2[1]);

	// Were coordinates found for each contig?
	if(coords1[0] == (int)path.getNumNodes() || coords2[0] == (int)path.getNumNodes())
	{
		start = path.getNumNodes();
		end = path.getNumNodes();
		// one cood missed
		return false;
	}
	
	size_t bestI = 0;
	size_t bestJ = 0;
	int best = path.getNumNodes();
	for(size_t i = 0; i <= 1; ++i)
		for(size_t j = 0; j <= 1; ++j)
		{
			int dist = abs(coords1[i] - coords2[j]);
			if(dist < best)
			{
				best = dist;
				bestI = i;
				bestJ = j;
			}
		}

	if(coords1[bestI] < coords2[bestJ])
	{
		start = coords1[bestI];
		end = coords2[bestJ];
	}
	else
	{
		start = coords2[bestJ];
		end = coords1[bestI];
	}
	
	return true;
	*/
	
}

static string toString(const ContigPath& path)
{
	size_t numNodes = path.getNumNodes();
	assert(numNodes > 0);
	MergeNode root = path.getNode(0);
	ostringstream s;
	s << root.id << (root.isRC ? '-' : '+');
	for (size_t i = 1; i < numNodes; ++i) {
		MergeNode mn = path.getNode(i);
		s << ',' << mn.id << (mn.isRC ? '-' : '+');
	}
	return s.str();
}

void mergePath(LinearNumKey cID, const ContigVec& sourceContigs,
		const ContigPath& currPath, int count, int kmer,
		FastaWriter* writer)
{
	if(gDebugPrint) std::cout << "Attempting to merge " << cID << "\n";
	if(gDebugPrint) std::cout << "Canonical path is: " << currPath << "\n"; 	
	string comment = toString(currPath);
	if (opt::verbose > 0)
		cout << comment << '\n';

	size_t numNodes = currPath.getNumNodes();

	MergeNode firstNode = currPath.getNode(0);
	Sequence merged = sourceContigs[firstNode.id].seq;
	if (firstNode.isRC)
		merged = reverseComplement(merged);
	assert(!merged.empty());

	for(size_t i = 1; i < numNodes; ++i) {
		MergeNode mn = currPath.getNode(i);
		if(gDebugPrint) std::cout << "	merging in " << mn.id << "(" << mn.isRC << ")\n";

		Sequence otherContig = sourceContigs[mn.id].seq;
		assert(!otherContig.empty());
		mergeSequences(merged, otherContig, (extDirection)0, mn.isRC, kmer);
	}

	writer->WriteSequence(merged, count, 0.0f, comment);
}


//
//
//
void mergeSequences(Sequence& rootContig, const Sequence& otherContig, extDirection dir, bool isReversed, size_t kmer)
{
	size_t overlap = kmer - 1;
	
	// should the slave be reversed?
	Sequence slaveSeq = otherContig;
	if(isReversed)
	{
		slaveSeq = reverseComplement(slaveSeq);
	}
	
	const Sequence* leftSeq;
	const Sequence* rightSeq;
	// Order the contigs
	if(dir == SENSE)
	{
		leftSeq = &rootContig;
		rightSeq = &slaveSeq;
	}
	else
	{
		leftSeq = &slaveSeq;
		rightSeq = &rootContig;
	}
	
	// Get the last k bases of the left and the first k bases of the right
	PackedSeq leftEnd = leftSeq->substr(leftSeq->length() - overlap, overlap);
	PackedSeq rightBegin = rightSeq->substr(0, overlap);
	
	// ensure that there is a legitimate k-1 overlap between these sequences	
	if(leftEnd != rightBegin)
	{
		printf("merge called data1: %s %s (%d, %d)\n", rootContig.c_str(), otherContig.c_str(), dir, isReversed);	
		printf("left end %s, right begin %s\n", leftEnd.decode().c_str(), rightBegin.decode().c_str());
		assert(leftEnd == rightBegin);
	}
	
	// TODO: make this in-place?
	// generate the merged sequence
	Sequence merged = *leftSeq;
	merged.append(rightSeq->substr(overlap));
	rootContig = merged;
}

void addPathNodesToList(MergeNodeList& list, ContigPath& path)
{
	size_t numNodes = path.getNumNodes();
	for(size_t idx = 0; idx < numNodes; idx++)
	{
		list.push_back(path.getNode(idx));
	}	
}

void parsePathLine(std::string pathLine, LinearNumKey& id, extDirection& dir, ContigPath& path)
{
	std::string discard;
	
	// discard the seperator
	std::stringstream pStream(pathLine);	
	pStream >> discard;
	
	// read in the root info
	std::string rootInfo;
	pStream >> rootInfo;

	int dirNum;
	readIDIntPair(rootInfo, id, dirNum);
	
	dir = (extDirection)dirNum;
	
	// discard the seperator
	pStream >> discard;
	
	pStream >> path;
}


void readIDIntPair(std::string str, LinearNumKey& id, int& i)
{
	std::stringstream ss(str);
	
	// read in the id
	std::string idLine;
	getline(ss, idLine, ',');
	
	id = convertContigIDToLinearNumKey(idLine);
	
	std::string intStr;
	ss >> intStr;
	
	std::stringstream convertor(intStr);
	convertor >> i;
}

