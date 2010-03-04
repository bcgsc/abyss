#include "config.h"
#include "Common/Options.h"
#include "AffixIterator.h"
#include "ContigPath.h"
#include "FastaReader.h"
#include "PairUtils.h"
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
#include <functional>
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
typedef pair<size_t, PathConsistencyStats> PathAlignment;

void readPathsFromFile(string pathFile, ContigPathMap& contigPathMap);
ContigPath* linkPaths(LinearNumKey id, ContigPathMap& paths,
		bool deleteSubsumed);
void mergePath(LinearNumKey cID, const ContigVec& sourceContigs,
		const ContigPath& mergeRecord, int count, int kmer,
		ostream& out);
void mergeSequences(Sequence& rootContig, const Sequence& otherContig, extDirection dir, bool isReversed, size_t kmer);
static PathAlignment align(const ContigPath& path1, ContigPath& path2,
		const ContigNode& pivot);
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
					affix_ostream_iterator<LinearNumKey>(cout, " "));
		else
			cout << " none";
		cout << '\n';
	}
	
	ostringstream ss;
	unsigned removed = 0;
	for (set<LinearNumKey>::const_iterator it = repeats.begin();
			it != repeats.end(); ++it) {
		if (paths.erase(*it) > 0) {
			ss << ' ' << idToString(*it);
			removed++;
		}
	}
	if (opt::verbose > 0 && removed > 0)
		cout << "Removing paths in repeats:" << ss.str() << '\n';
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
		FastaReader in(contigFile,
				FastaReader::KEEP_N | FastaReader::NO_FOLD_CASE);
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
			// Check if the current node to merge has any paths to/from it
			ContigPathMap::iterator findIter = paths.find(iter->id());
			if (findIter != paths.end()) {
				if (gDebugPrint)
					cout << ' ' << *iter << '\n';
				// Make the full path of the child node
				ContigPath childCanonPath = *findIter->second;

				if(gDebugPrint) cout << " ref: " << *refCanonical << '\n';
				if(gDebugPrint) cout << "  in: " << childCanonPath << '\n';

				PathAlignment a = align(*refCanonical, childCanonPath,
					*iter);
				bool validMerge = a.first > 0;
				size_t s2 = a.second.startP2, e2 = a.second.endP2;
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

template <class iterator>
static void skipAmbiguous(iterator& it1, iterator last1,
		iterator& it2, iterator last2)
{
	assert(it1 != last1);
	assert(it1->ambiguous());
	++it1;
	assert(it1 != last1);
	assert(!it1->ambiguous());

	assert(it2 != last2);
	it2 = find(it2, last2, *it1);
	if (it2 == last2)
		it1 = last1;
}

/** Find an equivalent region of the two specified paths, starting the
 * alignment at pos1 of path1 and pos2 of path2.
 * @return the alignment
 */
static PathAlignment align(
		const ContigPath& path1, const ContigPath& path2,
		ContigPath::const_iterator pivot1,
		ContigPath::const_iterator pivot2)
{
	ContigPath rpath2;
	bool flipped = pivot1->sense() != pivot2->sense();
	if (flipped) {
		rpath2 = path2;
		rpath2.reverseComplement();
		pivot2 = rpath2.begin() + (path2.end()-1 - pivot2);
	}
	const ContigPath& p1 = path1;
	const ContigPath& p2 = flipped ? rpath2 : path2;

	ContigPath::const_reverse_iterator rit1, rit2;
	for (rit1 = ContigPath::const_reverse_iterator(pivot1+1),
			rit2 = ContigPath::const_reverse_iterator(pivot2+1);
			rit1 != p1.rend() && rit2 != p2.rend();
			++rit1, ++rit2) {
		if (rit1->ambiguous() || rit2->ambiguous()) {
			if (rit1->ambiguous())
				skipAmbiguous<ContigPath::const_reverse_iterator>(
							rit1, p1.rend(), rit2, p2.rend());
			else
				skipAmbiguous<ContigPath::const_reverse_iterator>(
							rit2, p2.rend(), rit1, p1.rend());
			if (rit1 == p1.rend() || rit2 == p2.rend())
				break;
		}
		if (*rit1 != *rit2)
			break;
	}

	ContigPath::const_iterator it1, it2;
	for (it1 = pivot1, it2 = pivot2;
			it1 != p1.end() && it2 != p2.end();
			++it1, ++it2) {
		if (it1->ambiguous() || it2->ambiguous()) {
			if (it1->ambiguous())
				skipAmbiguous<ContigPath::const_iterator>(
						it1, p1.end(), it2, p2.end());
			else
				skipAmbiguous<ContigPath::const_iterator>(
						it2, p2.end(), it1, p1.end());
			if (it1 == p1.end() || it2 == p2.end())
				break;
		}
		if (*it1 != *it2)
			break;
	}

	if ((rit1 == p1.rend() || rit2 == p2.rend())
			&& (it1 == p1.end() || it2 == p2.end())) {
		// Found an alignment.
		assert((rit1 == p1.rend() && it2 == p2.end())
				|| (it1 == p1.end() && rit2 == p2.rend())
				|| (rit1 == p1.rend() && it1 == p1.end())
				|| (rit2 == p2.rend() && it2 == p2.end()));
		PathConsistencyStats a;
		a.startP1 = p1.rend() - rit1;
		a.endP1 = it1 - p1.begin() - 1;
		a.startP2 = p2.rend() - rit2;
		a.endP2 = it2 - p2.begin() - 1;
		a.flipped = flipped;
		a.duplicateSize = false;
		size_t count = max(a.endP1 - a.startP1 + 1,
				a.endP2 - a.startP2 + 1);
		return make_pair(count, a);
	} else
		return PathAlignment();
}

/** Return true if the IDs of the two ContigNodes are equal. */
static bool equalID(ContigNode a, ContigNode b)
{
	return a.id() == b.id();
}

/** Find an equivalent region of the two specified paths.
 * @param path2 is oriented to agree with path1
 * @return the alignment
 */
static PathAlignment align(const ContigPath& path1, ContigPath& path2,
		const ContigNode& pivot)
{
	assert(find(path1.begin(), path1.end(), pivot) != path1.end());
	ContigPath::iterator it2 = find(path2.begin(), path2.end(),
			pivot);
	if (it2 == path2.end())
		it2 = find(path2.begin(), path2.end(), ~pivot);
	assert(it2 != path2.end());
	assert(count(it2+1, path2.end(), pivot) == 0);

	map<size_t, PathConsistencyStats> pathAlignments;
	for (ContigPath::const_iterator it1 = find_if(
				path1.begin(), path1.end(),
				bind2nd(ptr_fun(equalID), pivot));
			it1 != path1.end();
			it1 = find_if(it1+1, path1.end(),
				bind2nd(ptr_fun(equalID), pivot))) {
		PathAlignment alignment = align(path1, path2, it1, it2);
		if (alignment.first == 0)
			continue;
		bool inserted = pathAlignments.insert(alignment).second;
		if (!inserted)
			pathAlignments[alignment.first].duplicateSize = true;
	}

	if (pathAlignments.empty()) {
		if (gDebugPrint)
			cout << " invalid " << pivot << '\n';
		return PathAlignment();
	}

	const PathAlignment& a = *pathAlignments.rbegin();

	// Sanity assert, at this point one of the low coordniates should
	// be zero and one of the high coordinates should be (size -1).
	size_t max1 = path1.size() - 1;
	size_t max2 = path2.size() - 1;
	assert(a.second.startP1 == 0 || a.second.startP2 == 0);
	assert(a.second.endP1 == max1 || a.second.endP2 == max2);

	// If either path aligns to the front and back of the other, it is
	// not a valid path.
	if (a.second.duplicateSize && a.first != min(max1+1, max2+1)) {
		if (gDebugPrint)
			cout << "Duplicate path match found\n";
		return PathAlignment();
	}

	if (a.second.flipped)
		path2.reverseComplement();

	assert(path1[a.second.startP1] == path2[a.second.startP2]
			|| (a.second.startP1 == 0 && a.second.startP2 == 0));
	assert(path1[a.second.endP1] == path2[a.second.endP2]
			|| (a.second.endP1 == max1 && a.second.endP2 == max2));
	return a;
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
