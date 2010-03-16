#include "config.h"
#include "Common/Options.h"
#include "AffixIterator.h"
#include "ContigPath.h"
#include "PairUtils.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <cstring> // for strerror
#include <deque>
#include <fstream>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
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
"Usage: " PROGRAM " [OPTION]... LEN PATH\n"
"Merge sequences of contigs IDs.\n"
"  LEN   lengths of the contigs\n"
"  PATH  sequences of contig IDs\n"
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

typedef map<LinearNumKey, ContigPath*> ContigPathMap;

/** Lengths of contigs. */
static vector<unsigned> g_contigLengths;

/** Return the length of this contig. */
unsigned ContigNode::length() const
{
	return ambiguous() ? m_id : g_contigLengths.at(id());
}

static ContigPath align(const ContigPath& p1, const ContigPath& p2,
		const ContigNode& pivot);

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
			transform(repeats.begin(), repeats.end(),
					affix_ostream_iterator<string>(cout, " "),
					idToString);
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

template <typename T> static const T& deref(const T* x)
{
	return *x;
}

/** Merge the specified path.
 * @param deleteSubsumed when true, remove paths that are entirely
 * subsumed by this path
 * Return a pointer to the merged path.
 */
static ContigPath* linkPaths(LinearNumKey id, ContigPathMap& paths,
		bool deleteSubsumed)
{
	ContigPath* path = deleteSubsumed ? paths[id]
		: new ContigPath(*paths[id]);
	if (gDebugPrint) {
		if (!deleteSubsumed)
			cout << "\n* " << ContigNode(id, false) << '\n';
		else
			cout << '\n' << ContigNode(id, false);
		cout << '\t' << *path << '\n';
	}

	set<ContigNode> seen;
	seen.insert(ContigNode(id, false));
	deque<ContigNode> mergeQ;
	for (ContigPath::const_iterator it = path->begin();
			it != path->end(); ++it)
		if (seen.insert(*it).second)
			mergeQ.push_back(*it);

	for (ContigNode pivot; !mergeQ.empty(); mergeQ.pop_front()) {
		pivot = mergeQ.front();
		ContigPathMap::iterator findIter = paths.find(pivot.id());
		if (findIter == paths.end())
			continue;
		ContigPath path2 = *findIter->second;
		if (pivot.sense())
			path2.reverseComplement();
		if (gDebugPrint)
			cout << pivot << '\t' << path2 << '\n';
		ContigPath consensus = align(*path, path2, pivot);
		if (consensus.empty())
			continue;
		if (deleteSubsumed) {
			/* If additional merges could be made at this point,
			 * something is wrong. We may need to delete all merged
			 * paths that exist for these paths and print the
			 * originals, but for now we keep both and print a
			 * warning.
			 */
			if (consensus != *path) {
				set<LinearNumKey> refKeys, childKeys;
				for (ContigPath::const_iterator it = path->begin();
						it != path->end(); it++)
					refKeys.insert(it->id());
				for (ContigPath::const_iterator it = path2.begin();
						it != path2.end(); it++)
					childKeys.insert(it->id());
				bool refIncludesChild
					= includes(refKeys.begin(), refKeys.end(),
							childKeys.begin(), childKeys.end());
				bool childIncludesRef
					= includes(childKeys.begin(), childKeys.end(),
							refKeys.begin(), refKeys.end());

				if (!refIncludesChild && !childIncludesRef) {
					/* A merge is now possible that wasn't made in the
					 * first pass. Ideally, this shouldn't happen, but
					 * can when a path is not extended due to its
					 * being an inverted repeat. Perhaps this
					 * restriction should be relaxed.
					 */
					if (gDebugPrint)
						cout << "\tnot merged\n";
				} else if (refIncludesChild && !childIncludesRef ) {
					if (gDebugPrint)
						cout << "\tremoved circular path\n";
					delete findIter->second;
					findIter->second = NULL;
					paths.erase(findIter);
				} else if (gDebugPrint)
					cout << "\twarning: possible circular path\n";
			} else {
				delete findIter->second;
				findIter->second = NULL;
				paths.erase(findIter);
			}
		} else {
			for (ContigPath::const_iterator it = path2.begin();
					it != path2.end(); ++it)
				if (seen.insert(*it).second)
					mergeQ.push_back(*it);

			path->swap(consensus);
			if (gDebugPrint)
				cout << '\t' << *path << '\n';
		}
	}
	return path;
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

static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

/** Read contig lengths. */
static vector<unsigned> readContigLengths(const string& path)
{
	vector<unsigned> lengths;
	ifstream in(path.c_str());
	assert_open(in, path);

	assert(g_contigIDs.empty());
	string id;
	unsigned len;
	while (in >> id >> len) {
		in.ignore(numeric_limits<streamsize>::max(), '\n');
		(void)g_contigIDs.serial(id);
		lengths.push_back(len);
	}
	assert(in.eof());
	assert(!lengths.empty());
	return lengths;
}

/** Read a set of paths from the specified file. */
static ContigPathMap readPaths(const string& filePath)
{
	ifstream in(filePath.c_str());
	assert_open(in, filePath);

	ContigPathMap paths;
	string idString;
	ContigPath path;
	while (in >> idString >> path) {
		LinearNumKey id = g_contigIDs.serial(idString);
		bool inserted = paths.insert(
				make_pair(id, new ContigPath(path))).second;
		assert(inserted);
		(void)inserted;
	}
	assert(in.eof());
	return paths;
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

	g_contigLengths = readContigLengths(argv[optind++]);
	g_contigIDs.lock();

	ContigPathMap originalPathMap = readPaths(argv[optind++]);

	removeRepeats(originalPathMap);

	ContigPathMap resultsPathMap;
	for (ContigPathMap::const_iterator iter = originalPathMap.begin();
			iter != originalPathMap.end(); ++iter)
		extendPaths(iter->first, originalPathMap, resultsPathMap);
	if (gDebugPrint)
		cout << '\n';

	removeRepeats(resultsPathMap);

	if (gDebugPrint)
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

/** Add the number of k-mer in two contigs. */
static unsigned addLength(unsigned addend, const ContigNode& contig)
{
	return addend + contig.length() - opt::k + 1;
}

/** Align the ambiguous region [it1, last1) to [it2, last2).
 * The end of the alignment is returned in it1 and it2 if there is
 * exactly one alignment. If there is more than one possible
 * alignment, the alignment is returned in it1 and the returned vector
 * of iterators.
 */
template <class iterator, class oiterator>
static vector<iterator> skipAmbiguous(iterator& it1, iterator last1,
		iterator& it2, iterator last2,
		oiterator out)
{
	(void)last1;
	assert(it1 != last1);
	assert(it1->ambiguous());
	assert(it1 + 1 != last1);
	assert(it2 != last2);

	unsigned ambiguous = it1->length();
	ContigNode needle = *++it1;
	assert(!needle.ambiguous());
	iterator it = find(it2, last2, needle);
	unsigned nmatches = count(it, last2, needle);

	vector<iterator> matches;
	switch (nmatches) {
	  case 0: {
		copy(it2, last2, out);
		unsigned unambiguous = accumulate(it2, last2, 0, addLength);
		it2 = last2;
		assert(ambiguous >= unambiguous);
		if (ambiguous > unambiguous)
			*out++ = ContigNode(ambiguous - unambiguous);
		break;
	  }
	  case 1:
		copy(it2, it, out);
		it2 = it;
		break;
	  default:
		matches.reserve(nmatches);
		for (it = find_if(it, last2,
					bind2nd(equal_to<ContigNode>(), needle));
				it != last2;
				it = find_if(it+1, last2,
					bind2nd(equal_to<ContigNode>(), needle)))
			matches.push_back(it);
		assert(matches.size() == nmatches);
	}
	return matches;
}

template <class iterator, class oiterator>
static vector<iterator> alignAmbiguous(iterator& it1, iterator last1,
		iterator& it2, iterator last2, int& which,
		oiterator out)
{
	which = it1->ambiguous() && it2->ambiguous()
		? (it1->length() > it2->length() ? 0 : 1)
		: it1->ambiguous() ? 0
		: it2->ambiguous() ? 1
		: -1;
	return which == -1 ? vector<iterator>()
		: which == 0 ? skipAmbiguous(it1, last1, it2, last2, out)
		: skipAmbiguous(it2, last2, it1, last1, out);
}

template <class iterator, class oiterator>
static bool align(iterator it1, iterator last1,
		iterator it2, iterator last2,
		oiterator out)
{
	assert(it1 != last1);
	assert(it2 != last2);
	for (; it1 != last1 && it2 != last2; ++it1, ++it2) {
		int which;
		vector<iterator> its = alignAmbiguous(it1, last1, it2, last2,
				which, out);
		if (!its.empty()) {
			// More than one match. Recurse on each option.
			assert(which == 0 || which == 1);
			for (typename vector<iterator>::iterator
					it = its.begin(); it != its.end(); ++it) {
				ContigPath consensus;
				if (align(which == 0 ? it1 : *it, last1,
							which == 1 ? it2 : *it, last2,
							back_inserter(consensus))) {
					copy(which == 0 ? it2 : it1, *it, out);
					copy(consensus.begin(), consensus.end(), out);
					return true;
				}
			}
			return false;
		}
		if (it1 == last1 || it2 == last2)
			break;

		if (*it1 != *it2)
			return false;
		*out++ = *it1;
	}

	assert(it1 == last1 || it2 == last2);
	copy(it1, last1, out);
	copy(it2, last2, out);
	return true;
}

/** Find an equivalent region of the two specified paths, starting the
 * alignment at pivot1 of path1 and pivot2 of path2.
 * @return the consensus sequence
 */
static ContigPath align(const ContigPath& p1, const ContigPath& p2,
		ContigPath::const_iterator pivot1,
		ContigPath::const_iterator pivot2)
{
	ContigPath::const_reverse_iterator
		rit1 = ContigPath::const_reverse_iterator(pivot1+1),
		rit2 = ContigPath::const_reverse_iterator(pivot2+1);
	ContigPath alignmentr;
	bool alignedr = align(rit1, p1.rend(), rit2, p2.rend(),
			back_inserter(alignmentr));

	ContigPath alignmentf;
	bool alignedf = align(pivot1, p1.end(), pivot2, p2.end(),
			back_inserter(alignmentf));

	ContigPath consensus;
	if (alignedr && alignedf) {
		// Found an alignment.
		assert(!alignmentf.empty());
		assert(!alignmentr.empty());
		consensus.reserve(alignmentr.size()-1 + alignmentf.size());
		consensus.assign(alignmentr.rbegin(), alignmentr.rend()-1);
		consensus.insert(consensus.end(),
				alignmentf.begin(), alignmentf.end());
	}
	return consensus;
}

/** Find an equivalent region of the two specified paths.
 * @return the consensus sequence
 */
static ContigPath align(
		const ContigPath& path1, const ContigPath& path2,
		const ContigNode& pivot)
{
	assert(find(path1.begin(), path1.end(), pivot) != path1.end());
	ContigPath::const_iterator it2 = find(path2.begin(), path2.end(),
			pivot);
	assert(it2 != path2.end());
	assert(count(it2+1, path2.end(), pivot) == 0);

	ContigPath consensus;
	for (ContigPath::const_iterator it1 = find_if(
				path1.begin(), path1.end(),
				bind2nd(equal_to<ContigNode>(), pivot));
			it1 != path1.end();
			it1 = find_if(it1+1, path1.end(),
				bind2nd(equal_to<ContigNode>(), pivot))) {
		consensus = align(path1, path2, it1, it2);
		if (!consensus.empty())
			return consensus;
	}

	if (gDebugPrint)
		cout << "\tinvalid\n";
	return consensus;
}
