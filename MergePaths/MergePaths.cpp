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
#if _OPENMP
# include <omp.h>
#endif

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
"  -j, --threads=N       use N parallel threads\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	static unsigned k;
	static string out;
	static int threads;
}

static const char shortopts[] = "j:k:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "out",         required_argument, NULL, 'o' },
	{ "threads",     required_argument,	NULL, 'j' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

typedef map<LinearNumKey, ContigPath> ContigPathMap;

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
		const ContigPath& path = pathIt->second;
		map<LinearNumKey, unsigned> count;
		for (ContigPath::const_iterator it = path.begin();
				it != path.end(); ++it)
			if (!it->ambiguous())
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

/** Extend the specified path as long as is unambiguously possible and
 * add the result to the specified container.
 */
static void extendPaths(LinearNumKey id,
		const ContigPathMap& paths, ContigPathMap& out)
{
	ContigPathMap::const_iterator pathIt = paths.find(id);
	assert(pathIt != paths.end());
	pair<ContigPathMap::iterator, bool> inserted;
	#pragma omp critical(out)
	inserted = out.insert(*pathIt);
	assert(inserted.second);
	ContigPath& path = inserted.first->second;

	if (gDebugPrint)
		#pragma omp critical(cout)
		cout << "\n* " << ContigNode(id, false) << '\n'
			<< '\t' << path << '\n';

	set<ContigNode> seen;
	seen.insert(ContigNode(id, false));
	deque<ContigNode> mergeQ;
	for (ContigPath::const_iterator it = path.begin();
			it != path.end(); ++it)
		if (seen.insert(*it).second)
			mergeQ.push_back(*it);

	for (ContigNode pivot; !mergeQ.empty(); mergeQ.pop_front()) {
		pivot = mergeQ.front();
		ContigPathMap::const_iterator path2It
			= paths.find(pivot.id());
		if (path2It == paths.end())
			continue;
		ContigPath path2 = path2It->second;
		if (pivot.sense())
			path2.reverseComplement();
		ContigPath consensus = align(path, path2, pivot);
		if (consensus.empty())
			continue;

		for (ContigPath::const_iterator it = path2.begin();
				it != path2.end(); ++it)
			if (seen.insert(*it).second)
				mergeQ.push_back(*it);

		path.swap(consensus);
		if (gDebugPrint)
			#pragma omp critical(cout)
			cout << '\t' << path << '\n';
	}
}

/** Remove paths subsumed by the specified path. */
static void removeSubsumedPaths(LinearNumKey id,
		ContigPathMap& paths)
{
	const ContigPath& path = paths[id];
	if (gDebugPrint)
		cout << '\n' << ContigNode(id, false)
			<< '\t' << path << '\n';

	for (ContigPath::const_iterator it = path.begin();
			it != path.end(); ++it) {
		ContigNode pivot = *it;
		if (pivot.id() == id)
			continue;
		ContigPathMap::iterator path2It = paths.find(pivot.id());
		if (path2It == paths.end())
			continue;
		ContigPath path2 = path2It->second;
		if (pivot.sense())
			path2.reverseComplement();
		ContigPath consensus = align(path, path2, pivot);
		if (consensus.empty())
			continue;
		assert(consensus == path);
		paths.erase(path2It);
	}
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
				make_pair(id, path)).second;
		assert(inserted);
		(void)inserted;
	}
	assert(in.eof());
	return paths;
}

#if _OPENMP
/** Store it in out and increment it.
 * @return true if out != last
 */
template<class T>
bool atomicInc(T& it, T last, T& out)
{
	#pragma omp critical(atomicInc)
	out = it == last ? it : it++;
	return out != last;
}
#endif

int main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'j': arg >> opt::threads; break;
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

#if _OPENMP
	if (opt::threads > 0)
		omp_set_num_threads(opt::threads);
#endif

	g_contigLengths = readContigLengths(argv[optind++]);
	g_contigIDs.lock();

	ContigPathMap originalPathMap = readPaths(argv[optind++]);

	removeRepeats(originalPathMap);

	ContigPathMap resultsPathMap;
#if _OPENMP
	ContigPathMap::iterator sharedIt = originalPathMap.begin();
	#pragma omp parallel
	for (ContigPathMap::iterator it;
			atomicInc(sharedIt, originalPathMap.end(), it);)
		extendPaths(it->first, originalPathMap, resultsPathMap);
#else
	for (ContigPathMap::const_iterator iter = originalPathMap.begin();
			iter != originalPathMap.end(); ++iter)
		extendPaths(iter->first, originalPathMap, resultsPathMap);
#endif
	originalPathMap.clear();
	if (gDebugPrint)
		cout << '\n';

	removeRepeats(resultsPathMap);

	if (gDebugPrint)
		cout << "\nRemoving redundant contigs\n";

	for (ContigPathMap::const_iterator iter = resultsPathMap.begin();
			iter != resultsPathMap.end(); ++iter)
		removeSubsumedPaths(iter->first, resultsPathMap);

	vector<ContigPath> uniquePaths;
	uniquePaths.reserve(resultsPathMap.size());
	for (ContigPathMap::const_iterator it = resultsPathMap.begin();
			it != resultsPathMap.end(); ++it)
		uniquePaths.push_back(it->second);
	sort(uniquePaths.begin(), uniquePaths.end());

	int id;
	istringstream ss(g_contigIDs.key(g_contigLengths.size() - 1));
	ss >> id;
	id++;

	ofstream fout(opt::out.c_str());
	ostream& out = opt::out.empty() ? cout : fout;
	assert(out.good());
	for (vector<ContigPath>::const_iterator it
				= uniquePaths.begin();
			it != uniquePaths.end(); ++it)
		out << id++ << '\t' << *it << '\n';
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
	assert(it1 != last1);
	assert(it1->ambiguous());
	assert(it1 + 1 != last1);
	assert(it2 != last2);

	// Find a seed for the alignment.
	unsigned ambiguous1 = it1->length();
	iterator it1e = ++it1;
	for (iterator it = it1; it != last1; ++it) {
		if (it->ambiguous())
			continue;
		iterator it2e = find(it2, last2, *it);
		if (it2e != last2) {
			it1e = it;
			break;
		}
	}
	assert(!it1e->ambiguous());
	iterator it2e = find(it2, last2, *it1e);
	unsigned nmatches = count(it2e, last2, *it1e);

	vector<iterator> matches;
	switch (nmatches) {
	  case 0: {
		// Unable to find the seed in path2. Check whether the
		// remainder of path2 fits entirely within the gap of path1.
		unsigned unambiguous2 = accumulate(it2, last2, 0, addLength);
		if (ambiguous1 < unambiguous2) {
			// The size of the seqeuence in path2 is larger than the
			// gap in path1. No alignment.
			return matches;
		}
		copy(it2, last2, out);
		it1 = it1e;
		it2 = last2;
		if (ambiguous1 > unambiguous2)
			*out++ = ContigNode(ambiguous1 - unambiguous2);
		break;
	  }
	  case 1:
		// The seed occurs exactly once in path2.
		if (it1 == it1e) {
			// path2 completely fills the gap in path1.
			copy(it2, it2e, out);
			it2 = it2e;
		} else {
			// The gaps of path1 and path2 overlap.
			iterator it2a = it2e - 1;
			if (it2e == it2 || !it2a->ambiguous()) {
				// The two paths do not agree. No alignment.
				return matches;
			}

			unsigned ambiguous2 = it2a->length();
			unsigned unambiguous1 = accumulate(it1, it1e,
					0, addLength);
			assert(ambiguous2 >= unambiguous1);

			unsigned unambiguous2 = accumulate(it2, it2a,
					0, addLength);
			assert(ambiguous1 >= unambiguous2);

			unsigned n = min(ambiguous2 - unambiguous1,
					ambiguous1 - unambiguous2);

			out = copy(it2, it2a, out);
			if (n > 0)
				*out++ = ContigNode(n);
			out = copy(it1, it1e, out);
			it1 = it1e;
			it2 = it2e;
		}
		break;
	  default:
		// The seed occurs more than once in path2. Return all the
		// matches so that our caller may iterate over them.
		assert(it1 == it1e);
		matches.reserve(nmatches);
		for (iterator it = find_if(it2e, last2,
					bind2nd(equal_to<ContigNode>(), *it1e));
				it != last2;
				it = find_if(it+1, last2,
					bind2nd(equal_to<ContigNode>(), *it1e)))
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
	if (gDebugPrint)
		#pragma omp critical(cout)
		cout << pivot << '\t' << path2 << '\n';

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
		#pragma omp critical(cout)
		cout << "\tinvalid\n";
	return consensus;
}
