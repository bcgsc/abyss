#include "config.h"
#include "Common/Options.h"
#include "ContigPath.h"
#include "FastaReader.h"
#include "PairUtils.h"
#include "PrefixIterator.h"
#include "Sense.h"
#include "Sequence.h"
#include "StringUtil.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>

using namespace std;

#define PROGRAM "MergePaths"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Jared Simpson and Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... [CONTIG] PATH\n"
"Merge paths and contigs. If CONTIG is specified, the output is\n"
"FASTA and merged paths otherwise.\n"
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
	static string out;
}

static const char shortopts[] = "k:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "out",         required_argument, NULL, 'o' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

struct PathConsistencyStats {
	size_t startP1;
	size_t endP1;
	size_t startP2;
	size_t endP2;
	bool flipped;
	bool duplicateSize;
};

typedef list<MergeNode> MergeNodeList;
typedef map<LinearNumKey, ContigPath*> ContigPathMap;

struct Contig {
	string id;
	Sequence seq;
	unsigned coverage;

	Contig(const string& id, const Sequence& seq, unsigned coverage)
		: id(id), seq(seq), coverage(coverage) { }

	friend ostream& operator <<(ostream& out, const Contig& o)
	{
		return out << '>' << o.id << ' '
			<< o.seq.length() << ' ' << o.coverage << '\n'
			<< o.seq << '\n';
	}
};

typedef vector<Contig> ContigVec;

void readPathsFromFile(string pathFile, ContigPathMap& contigPathMap);
ContigPath* linkPaths(LinearNumKey id, ContigPathMap& paths,
		bool deleteSubsumed);
void mergePath(LinearNumKey cID, const ContigVec& sourceContigs,
		const ContigPath& mergeRecord, int count, int kmer,
		ostream& out);
void mergeSequences(Sequence& rootContig, const Sequence& otherContig, extDirection dir, bool isReversed, size_t kmer);
bool extractMinCoordSet(LinearNumKey anchor, ContigPath& path, vector<size_t>& coords);
bool checkPathConsistency(LinearNumKey path1Root, LinearNumKey path2Root, ContigPath& path1, ContigPath& path2, size_t& startP1, size_t& endP1, size_t& startP2, size_t& endP2);
void addPathNodesToList(MergeNodeList& list, ContigPath& path);

static bool gDebugPrint;

/** Return all contigs that are tandem repeats, identified as those
 * contigs that appear more than once in a single path.
 */
static set<LinearNumKey> findRepeats(const ContigPathMap& paths)
{
	set<LinearNumKey> repeats;
	for (ContigPathMap::const_iterator pathIt = paths.begin();
			pathIt != paths.end(); ++pathIt) {
		const ContigPath& path = *pathIt->second;
		map<LinearNumKey, unsigned> count;
		for (ContigPath::const_iterator it = path.begin();
				it != path.end(); ++it)
			count[it->id()]++;
		for (map<LinearNumKey, unsigned>::const_iterator
				it = count.begin(); it != count.end(); ++it)
			if (it->second > 1)
				repeats.insert(it->first);
	}
	return repeats;
}

/** Convert a numeric contig ID to a string. */
static const string& idToString(unsigned id)
{
	return g_contigIDs.key(id);
}

/** Remove tandem repeats from the set of paths. */
static void removeRepeats(ContigPathMap& paths)
{
	set<LinearNumKey> repeats = findRepeats(paths);
	if (gDebugPrint) {
		cout << "Repeats:";
		if (!repeats.empty())
			copy(repeats.begin(), repeats.end(),
					prefix_ostream_iterator<LinearNumKey>(cout, " "));
		else
			cout << " none";
		cout << '\n';
	}
	cout << "Removing paths in repeats:";
	unsigned removed = 0;
	for (set<LinearNumKey>::const_iterator it = repeats.begin();
			it != repeats.end(); ++it) {
		if (paths.erase(*it) > 0) {
			cout << ' ' << idToString(*it);
			removed++;
		}
	}
	if (removed == 0)
		cout << " none";
	cout << '\n';
}

static set<size_t> getContigIDs(const vector<ContigPath>& paths)
{
	set<size_t> seen;
	for (vector<ContigPath>::const_iterator it = paths.begin();
			it != paths.end(); it++) {
		size_t nodes = it->size();
		for (size_t i = 0; i < nodes; i++)
			seen.insert((*it)[i].id());
	}
	return seen;
}

template <typename T> static const T& deref(const T* x)
{
	return *x;
}

/** Extend the specified path as long as is unambiguously possible and
 * add the result to the specified container. */
static void extendPaths(LinearNumKey id,
		ContigPathMap& paths, ContigPathMap& out)
{
	bool inserted = out.insert(make_pair(id,
				linkPaths(id, paths, false))).second;
	assert(inserted);
	(void)inserted;
}

/** Remove paths subsumed by the specified path. */
static void removeSubsumedPaths(LinearNumKey id,
		ContigPathMap& paths)
{
	linkPaths(id, paths, true);
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

	if (argc - optind > 1) {
		if (opt::k <= 0) {
			cerr << PROGRAM ": missing -k,--kmer option\n";
			die = true;
		}

		if (opt::out.empty()) {
			cerr << PROGRAM ": " << "missing -o,--out option\n";
			die = true;
		}
	}

	if (argc - optind < 1) {
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

	const char* contigFile = argc - optind == 1 ? NULL
		: argv[optind++];
	string pathFile(argv[optind++]);

	ContigVec contigVec;
	if (contigFile != NULL) {
		FastaReader in(contigFile, FastaReader::KEEP_N);
		for (FastaRecord rec; in >> rec;) {
			istringstream ss(rec.comment);
			unsigned length, coverage = 0;
			ss >> length >> coverage;
			LinearNumKey id = g_contigIDs.serial(rec.id);
			assert(id == contigVec.size());
			(void)id;
			contigVec.push_back(Contig(rec.id, rec.seq, coverage));
		}
		g_contigIDs.lock();
		assert(in.eof());
		assert(!contigVec.empty());
		opt::colourSpace = isdigit(contigVec[0].seq[0]);
	}

	// Read the paths file
	ContigPathMap originalPathMap, resultsPathMap;
	readPathsFromFile(pathFile, originalPathMap);

	removeRepeats(originalPathMap);

	// link the paths together
	for (ContigPathMap::const_iterator iter = originalPathMap.begin();
			iter != originalPathMap.end(); ++iter) {
		extendPaths(iter->first, originalPathMap, resultsPathMap);
		if (gDebugPrint)
			cout << "Merged path: "
				<< *resultsPathMap[iter->first] << '\n';
	}

	if (opt::verbose > 0)
		cout << "\nRemoving redundant contigs\n";

	for (ContigPathMap::const_iterator iter = resultsPathMap.begin();
			iter != resultsPathMap.end(); ++iter)
		removeSubsumedPaths(iter->first, resultsPathMap);

	set<ContigPath*> uniquePtr;
	for (ContigPathMap::const_iterator it = resultsPathMap.begin();
			it != resultsPathMap.end(); ++it)
		uniquePtr.insert(it->second);

	// Sort the set of unique paths by the path itself rather than by
	// pointer. This ensures that the order of the contig IDs does not
	// depend on arbitrary pointer values.
	vector<ContigPath> uniquePaths;
	uniquePaths.reserve(uniquePtr.size());
	transform(uniquePtr.begin(), uniquePtr.end(),
			back_inserter(uniquePaths), deref<ContigPath>);
	sort(uniquePaths.begin(), uniquePaths.end());

	if (contigVec.empty()) {
		ofstream fout(opt::out.c_str());
		ostream& out = opt::out.empty() ? cout : fout;
		assert(out.good());
		unsigned pathID = 0;
		for (vector<ContigPath>::const_iterator it
					= uniquePaths.begin();
				it != uniquePaths.end(); ++it)
			out << pathID++ << '\t' << *it << '\n';
		assert(out.good());
		return 0;
	}

	if (opt::verbose > 0)
		cout << "\nMerging contigs\n";

	ofstream out(opt::out.c_str());
	set<size_t> seen = getContigIDs(uniquePaths);
	float minCov = numeric_limits<float>::infinity(),
		  minCovUsed = numeric_limits<float>::infinity();
	for (size_t i = 0; i < contigVec.size(); i++) {
		const Contig& contig = contigVec[i];
		bool used = seen.count(i) > 0;
		if (!used)
			out << contig;
		if (contig.coverage > 0) {
			assert((int)contig.seq.length() - opt::k + 1 > 0);
			float cov = (float)contig.coverage
				/ (contig.seq.length() - opt::k + 1);
			minCov = min(minCov, cov);
			if (used)
				minCovUsed = min(minCovUsed, cov);
		}
	}

	cout << "The minimum coverage of single-end contigs is "
		<< minCov << ".\n"
		<< "The minimum coverage of merged contigs is "
		<< minCovUsed << ".\n";
	if (minCov < minCovUsed)
		cout << "Consider increasing the coverage threshold "
			"parameter, c, to " << minCovUsed << ".\n";

	stringstream s(g_contigIDs.key(contigVec.size() - 1));
	int id;
	s >> id;
	id++;
	for (vector<ContigPath>::const_iterator it = uniquePaths.begin();
			it != uniquePaths.end(); ++it)
		mergePath(it->front().id(), contigVec, *it, id++,
				opt::k, out);
	assert(out.good());

	return 0;
}

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

void readPathsFromFile(string pathFile, ContigPathMap& contigPathMap)
{
	ifstream pathStream(pathFile.c_str());
	assert_open(pathStream, pathFile);

	string line;
	while (getline(pathStream, line)) {
		istringstream s(line);
		string sID;
		ContigPath path;
		s >> sID >> path;
		assert(s.eof());

		LinearNumKey id = g_contigIDs.serial(sID);
		bool inserted = contigPathMap.insert(make_pair(id,
			new ContigPath(path))).second;
		assert(inserted);
		(void)inserted;
	}
	assert(pathStream.eof());

	pathStream.close();
}

/** Merge the specified path.
 * @param deleteSubsumed when true, remove paths that are entirely
 * subsumed by this path
 * Return a pointer to the merged path.
 */
ContigPath* linkPaths(LinearNumKey id, ContigPathMap& paths,
		bool deleteSubsumed)
{
	ContigPath* refCanonical = deleteSubsumed
		? paths[id] : new ContigPath(*paths[id]);
	if (gDebugPrint)
		cout << "\n* " << idToString(id) << ": "
			<< *refCanonical << '\n';

	// Build the initial list of nodes to attempt to merge in
	MergeNodeList mergeInList;
	addPathNodesToList(mergeInList, *refCanonical);

	MergeNodeList::iterator iter = mergeInList.begin();
	while(!mergeInList.empty()) {
		if(iter->id() != id) {
			if (gDebugPrint) cout << ' ' << *iter << '\n';

			// Check if the current node to merge has any paths to/from it
			ContigPathMap::iterator findIter = paths.find(iter->id());
			if (findIter != paths.end()) {
				// Make the full path of the child node
				ContigPath childCanonPath = *findIter->second;

				if(gDebugPrint) cout << " ref: " << *refCanonical << '\n';
				if(gDebugPrint) cout << "  in: " << childCanonPath << '\n';

				size_t s1, s2, e1, e2;
				bool validMerge = checkPathConsistency(id, iter->id(),
					*refCanonical, childCanonPath, s1, e1, s2, e2);

				if(validMerge && deleteSubsumed) {
					// If additional merges could be made at this
					// point, something is wrong. We may need to delete
					// all merged paths that exist for these paths and
					// print the originals, but for now we keep both
					// and print a warning.
					if (s2 != 0 || e2+1 != childCanonPath.size()) {
						set<LinearNumKey> refKeys, childKeys;
						for (ContigPath::const_iterator it = refCanonical->begin();
								it != refCanonical->end(); it++)
							refKeys.insert(it->id());
						for (ContigPath::const_iterator it = childCanonPath.begin();
								it != childCanonPath.end(); it++)
							childKeys.insert(it->id());
						bool refIncludesChild =
							includes(refKeys.begin(), refKeys.end(),
									childKeys.begin(), childKeys.end());
						bool childIncludesRef =
							includes(childKeys.begin(), childKeys.end(),
									refKeys.begin(), refKeys.end());

						if (!refIncludesChild && !childIncludesRef) {
							/* A merge is now possible that wasn't
							 * made in the first pass. Ideally, this
							 * shouldn't happen, but can when a path
							 * is not extended due to its being an
							 * inverted repeat. Perhaps this
							 * restriction should be relaxed.
							 */
							if (gDebugPrint)
								cout << " not merging: " << childCanonPath << '\n';
						} else if (refIncludesChild && !childIncludesRef ) {
							if(gDebugPrint)
								cout << " removing circular: " << childCanonPath << '\n';
							delete findIter->second;
							findIter->second = NULL;
							paths.erase(findIter);
						} else if (gDebugPrint)
							cout << " warning: possible circular paths\n";
					} else {
						if (gDebugPrint)
							cout << " del: " << childCanonPath << '\n';
						delete findIter->second;
						findIter->second = NULL;
						paths.erase(findIter);
					}
				} else if (validMerge) {
					// Extract the extra nodes from the child path that can be added in
					ContigPath prependNodes(&childCanonPath[0],
							&childCanonPath[s2]);
					ContigPath appendNodes(&childCanonPath[e2+1],
							&childCanonPath[childCanonPath.size()]);

					// Add the nodes to the list of contigs to try to merge in
					addPathNodesToList(mergeInList, prependNodes);
					addPathNodesToList(mergeInList, appendNodes);

					// Add the nodes to the ref contig
					refCanonical->insert(refCanonical->begin(),
							prependNodes.begin(), prependNodes.end());
					refCanonical->insert(refCanonical->end(),
							appendNodes.begin(), appendNodes.end());

					if(gDebugPrint) cout << " new: " << *refCanonical << '\n';
				}
			}
		}
		// Erase the iterator and move forward
		mergeInList.erase(iter++);
	}
	return refCanonical;
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

	assert(path1.size() != 0 && path2.size() != 1);

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
	size_t max1 = path1.size() - 1;
	size_t max2 = path2.size() - 1;
	map<size_t, PathConsistencyStats> pathAlignments;

	for (unsigned i = 0; i < coords1.size(); i++) {
		for (unsigned j = 0; j < coords2.size(); j++) {
			startP1 = coords1[i];
			endP1 = coords1[i];
			if (flipped) {
				startP2 = max2 - coords2[j];
				endP2 = max2 - coords2[j];
			} else {
				startP2 = coords2[j];
				endP2 = coords2[j];
			}

			if(path1[startP1].sense() != path2[startP2].sense()) {
				path2.reverseComplement();
				flipped = !flipped;
				startP2 = max2 - startP2;
				endP2 = max2 - endP2;
			}

			lowValid = true;
			while(1) {
				if(path1[startP1].id() != path2[startP2].id()) {
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
				if(path1[endP1].id() != path2[endP2].id()) {
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
			if (lowValid && highValid) {
				size_t count = endP1 - startP1;
				assert(endP2 - startP2 == count);
				if (pathAlignments.find(count) == pathAlignments.end()) {
					PathConsistencyStats& pathAlignment = pathAlignments[count];
					pathAlignment.startP1 = startP1;
					pathAlignment.endP1 = endP1;
					pathAlignment.startP2 = startP2;
					pathAlignment.endP2 = endP2;
					pathAlignment.flipped = flipped;
					pathAlignment.duplicateSize = false;
				} else
					pathAlignments[count].duplicateSize = true;
			}
		}
	}

	// Check if there was an actual mismatch in the nodes
	if(pathAlignments.empty()) {
		if(gDebugPrint) printf("Invalid path match!\n");
		if(gDebugPrint) cout << "Path1 (" << path1Root << ") " << path1 << '\n';
		if(gDebugPrint) cout << "Path2 (" << path2Root << ") " << path2 << '\n';
		return false;
	}

	//printf("Final coords: [%zu-%zu] [%zu-%zu]\n", startP1, endP1, startP2, endP2);

	map<size_t, PathConsistencyStats>::const_iterator biggestIt =
		pathAlignments.end();
	--biggestIt;

	// Sanity assert, at this point one of the low coordniates should be zero and one of the high coordinates should be (size -1)
	assert(biggestIt->second.startP1 == 0 || biggestIt->second.startP2 == 0);
	assert(biggestIt->second.endP1 == max1 || biggestIt->second.endP2 == max2);

	// If either path aligns to the front and back of the other, it is
	// not a valid path.
	size_t count = biggestIt->first;
	if (biggestIt->second.duplicateSize
			&& biggestIt->first != min(max1, max2)) {
		if (gDebugPrint) printf("Duplicate path match found\n");
		return false;
	}

	startP1 = biggestIt->second.startP1;
	endP1 = biggestIt->second.endP1;
	startP2 = biggestIt->second.startP2;
	endP2 = biggestIt->second.endP2;
	if (biggestIt->second.flipped != flipped)
		path2.reverseComplement();

	for(size_t c = 0; c < count; ++c) {
		if(path1[startP1 + c].id() != path2[startP2 + c].id()) {
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
	size_t maxIdx = path.size();
	for(size_t idx = 0; idx < maxIdx; ++idx) {
		size_t tIdx = maxIdx - idx - 1;
		if(path[tIdx].id() == anchor)
			coords.push_back(tIdx);
	}

	if(coords.empty()) // anchor coord not found
		return false;

	return true;

	/*
	printf("	found %zu %zu %zu %zu\n", coords1[0], coords1[1], coords2[0], coords2[1]);

	// Were coordinates found for each contig?
	if (coords1[0] == (int)path.size()
			|| coords2[0] == (int)path.size()) {
		start = path.size();
		end = path.size();
		// one cood missed
		return false;
	}

	size_t bestI = 0;
	size_t bestJ = 0;
	int best = path.size();
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

/** Return a string representation of the specified object. */
template<typename T> static string toString(T x)
{
	ostringstream s;
	s << x;
	return s.str();
}

/** Return a string representation of the specified path. */
static string toString(const ContigPath& path, char sep)
{
	assert(!path.empty());
	ostringstream s;
	ContigPath::const_iterator it = path.begin();
	s << *it++;
	for (; it != path.end(); ++it)
		s << sep << *it;
	return s.str();
}

void mergePath(LinearNumKey cID, const ContigVec& sourceContigs,
		const ContigPath& currPath, int count, int kmer,
		ostream& out)
{
	if (gDebugPrint)
		cout << "Merging " << idToString(cID) << ": "
			<< currPath << '\n';
	if (opt::verbose > 0)
		cout << currPath << '\n';

	size_t numNodes = currPath.size();
	MergeNode firstNode = currPath[0];
	const Contig& firstContig = sourceContigs[firstNode.id()];
	Sequence merged = firstContig.seq;
	unsigned coverage = firstContig.coverage;
	if (firstNode.sense())
		merged = reverseComplement(merged);
	assert(!merged.empty());

	for(size_t i = 1; i < numNodes; ++i) {
		MergeNode mn = currPath[i];
		const Contig& contig = sourceContigs[mn.id()];
		assert(!contig.seq.empty());
		mergeSequences(merged, contig.seq, (extDirection)0, mn.sense(),
				kmer);
		coverage += contig.coverage;
	}

	ostringstream s;
	s << merged.length() << ' ' << coverage << ' '
		<< toString(currPath, ',');
	out << FastaRecord(toString(count), s.str(), merged);
}

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

	Sequence leftEnd(leftSeq->substr(leftSeq->length() - overlap,
				overlap));
	Sequence rightBegin(rightSeq->substr(0, overlap));
	if (leftEnd != rightBegin) {
		printf("merge called data1: %s %s (%d, %d)\n",
				rootContig.c_str(), otherContig.c_str(),
				dir, isReversed);
		printf("left end %s, right begin %s\n",
				leftEnd.c_str(), rightBegin.c_str());
		assert(leftEnd == rightBegin);
	}

	rootContig = *leftSeq + rightSeq->substr(overlap);
}

void addPathNodesToList(MergeNodeList& list, ContigPath& path)
{
	size_t numNodes = path.size();
	for(size_t idx = 0; idx < numNodes; idx++)
		list.push_back(path[idx]);
}
