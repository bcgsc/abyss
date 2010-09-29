/**
 *
 * Resolve ambiguity ("N"s) in paths after MergePaths
 * by using pairwise/multiple sequence alignment
 * (pairwise: Needleman-Wunsch)
 * (multi: based on Dialign program by A.R.Subramanian)
 *
 * Author: Rong She (rshe@bcgsc.ca)
 *
 * Last Modified: July 2010
 */

#include "dialign.h"
#include "config.h"
#include "Common/Options.h"
#include "ConstString.h"
#include "ConstrainedSearch.h"
#include "ContigGraph.h"
#include "ContigGraphAlgorithms.h"
#include "ContigNode.h"
#include "ContigPath.h"
#include "Dictionary.h"
#include "FastaReader.h"
#include "StringUtil.h"
#include "needleman_wunsch.h"
#include "smith_waterman.h"
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cstring> // for strerror
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "PathConsensus"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Rong She and Shaun Jackman.\n"
"\n"
"Copyright 2010 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... FASTA ADJ PATH\n"
"Align sequences of ambiguous paths and output a consensus\n"
"sequence.\n"
"  FASTA  contigs in FASTA format\n"
"  ADJ    contig adjacency graph\n"
"  PATH   paths of these contigs\n"
"\n"
" Options:\n"
"  -k, --kmer=N          k-mer size\n"
"  -o, --out=FILE        output contig paths to FILE\n"
"  -s, --consensus=FILE  output consensus sequences to FILE\n"
"  -a, --branches=N      maximum number of sequences to align\n"
"                        default: 4\n"
"  -p, --identity=REAL   minimum identity, default: 0.9\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
" DIALIGN-TX options:\n"
"  -D, --dialign-d=N     dialign debug level, default: 0\n"
"  -M, --dialign-m=FILE  score matrix, default: dna_matrix.scr\n"
"  -P, --dialign-p=FILE  diagonal length probability distribution\n"
"                        default: dna_diag_prob_100_exp_550000\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties
	static string out;
	static string consensusPath;
	static float identity = 0.9;
	static unsigned numBranches = 4;
	static int dialign_debug;
	static string dialign_score;
	static string dialign_prob;

	/** The the number of bases to continue the constrained search of
	 * the graph beyond the size of the ambiguous gap in the path.
	 */
	static const unsigned allowedError = 6;
}

static const char shortopts[] = "k:o:s:a:p:vD:M:P:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "kmer",        required_argument, NULL, 'k' },
	{ "out",         required_argument, NULL, 'o' },
	{ "consensus",   required_argument, NULL, 's' },
	{ "branches",    required_argument, NULL, 'a' },
	{ "identity",    required_argument, NULL, 'p' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ "dialign-d",   required_argument, NULL, 'D' },
	{ "dialign-m",   required_argument, NULL, 'M' },
	{ "dialign-p",   required_argument, NULL, 'P' },
	{ NULL, 0, NULL, 0 }
};

static struct {
	unsigned numAmbPaths;
	unsigned numMerged;
	unsigned numNoSolutions;
	unsigned numTooManySolutions;
	unsigned tooComplex;
	unsigned notMerged;
} stats;

struct AmbPathConstraint {
	ContigNode	source;
	ContigNode	dest;
	int		dist;

	AmbPathConstraint(const ContigNode& src, const ContigNode& dst,
			int d)
		: source(src), dest(dst), dist(d) {}

	bool operator ==(const AmbPathConstraint& a) const
	{
		return source == a.source && dest == a.dest
			&& dist == a.dist;
	}

	bool operator <(const AmbPathConstraint& a) const
	{
		return source < a.source
			|| (source == a.source && dest < a.dest)
			|| (source == a.source && dest == a.dest
				&& dist < a.dist);
	}
};

typedef ContigPath Path;
typedef vector<Path> ContigPaths;
typedef map<AmbPathConstraint, ContigPath> AmbPath2Contig;

typedef vector<const_string> Contigs;
static Contigs g_contigs;
AmbPath2Contig g_ambpath_contig;
static set< pair<ContigNode, ContigNode> > g_edges_irregular;

#if 0
static double g_coverage_mean;
static double g_coverage_variance;
#endif

/* commonly-used utility function */
static void assert_open(ifstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
}

/** Return the sequence of the specified contig node. The sequence
 * may be ambiguous or reverse complemented.
 */
static const Sequence getSequence(ContigNode id)
{
	if (id.ambiguous()) {
		string s(id.ambiguousSequence());
		if (s.length() < opt::k)
			transform(s.begin(), s.end(), s.begin(), ::tolower);
		return string(opt::k - 1, 'N') + s;
	} else {
		string seq(g_contigs[id.id()]);
		return id.sense() ? reverseComplement(seq) : seq;
	}
}

/** Read contig paths from the specified file.
 * @param ids [out] the string ID of the paths
 * @param isAmb [out] whether the path contains a gap
 */
static ContigPaths readPath(const string& inPath,
	vector<string>& ids, vector<bool>& isAmb)
{
	assert(ids.empty()); //this seed is contigID of the path
	assert(isAmb.empty());
	assert(g_ambpath_contig.empty());
	ifstream fin(inPath.c_str());
	if (opt::verbose > 0)
		cerr << "Reading `" << inPath << "'..." << endl;
	if (inPath != "-")
		assert_open(fin, inPath);
	istream& in = inPath == "-" ? cin : fin;

	ContigPaths paths;
	string id;
	Path path;
	Path::iterator prev, next;
	while (in >> id >> path) {
		paths.push_back(path);
		ids.push_back(id);

		if (path.size() <= 2) {
			isAmb.push_back(false);
			continue;
		}

		/* contig with Ns must have prev and next */
		bool cur_is_amb = false;
		prev = path.begin();
		next = path.begin();
		Path::iterator it = path.begin();
		it++; next++; next++;
		for (; next != path.end(); it++, prev++, next++) {
			if (it->ambiguous()) {
				cur_is_amb = true;
				/* add the entry to g_ambpath_contig
				(value 0: unprocessed/no proper paths) */
				g_ambpath_contig.insert(AmbPath2Contig::value_type(
					AmbPathConstraint(*prev, *next, -it->id()),
					ContigPath()));
			}
		}
		isAmb.push_back(cur_is_amb);
	}
	assert(in.eof());
	return paths;
}

#if 0
/** Read ambiguous path info from SimpleGraph verbose output.
 * format:
 * 275+ 1181+ 1287- 614+
 * 275+ 1198+ 1287- 614+
 * 275+ 73N 1287- 614+
 *
 * Store path info in "g_amb_paths"
 */
static void readAmbPaths(const string& ambPath)
{
	ifstream fin(ambPath.c_str());
	if (opt::verbose > 0)
		cerr << "Reading `" << ambPath << "'..." << endl;
	if (ambPath != "-")
		assert_open(fin, ambPath);
	istream& in = ambPath == "-" ? cin : fin;

	vector<Path> paths;
	string id;
	Path path;
	bool isAmbPath;
	while (in >> path) {
		isAmbPath = false;
		for (Path::iterator it = path.begin();
			it != path.end(); it++)
			if (it->ambiguous()) {
				isAmbPath = true;
				assert(!paths.empty());
				g_amb_paths.insert(
					map<Path, vector<Path> >::value_type(
					path, paths));
				paths.clear();
				break;
			}
		if (!isAmbPath)
			paths.push_back(path);
	}
	assert(in.eof());
}
#endif

/** Mark every contig in path as seen. */
static void markSeen(vector<bool>& seen, const ContigPath& path,
		bool flag)
{
	for (Path::const_iterator it = path.begin();
			it != path.end(); ++it)
		if (!it->ambiguous() && it->id() < seen.size())
			seen[it->id()] = flag;
}

/** Mark every contig in paths as seen. */
static void markSeen(vector<bool>& seen, const vector<Path>& paths,
		bool flag)
{
	for (vector<Path>::const_iterator it = paths.begin();
			it != paths.end(); ++it)
		markSeen(seen, *it, flag);
}

#if 0
/* read in adj file, store which contigs have irregular overlaps
 * (>=k-1 or non-exact match)
 */
static void ReadAdj(const string& irregularAdj)
{
	ifstream in_adj(irregularAdj.c_str());
	assert_open(in_adj, irregularAdj);
	while (!in_adj.eof()) {
		char line[65536];
		in_adj.getline(line, 65536);
		string s = line;
		vector<string> strs;
		size_t start_pos = 0, delimit_pos;
		while (start_pos < s.length()
			&& (delimit_pos = s.find(';', start_pos))
			!= string::npos) {
			//get first, second parts
			if (delimit_pos > start_pos + 1)
				strs.push_back(s.substr(start_pos,
					delimit_pos-start_pos));
			else
				strs.push_back("");
			start_pos = delimit_pos+1;
		}
		//ignore lines that do not have two ';'s
		if (strs.size() < 2)
			continue;
		//get the last part
		if (start_pos + 1 < s.length())
			strs.push_back(s.substr(start_pos));
		else
			strs.push_back("");
		assert(strs.size() == 3);
		size_t space_pos;
		assert ((space_pos = strs[0].find(' '))
			!= string::npos);
		string t_id = strs[0].substr(0,space_pos);
		if (!strs[1].empty()) {
			start_pos = 0;
			while (start_pos < strs[1].length()
				&& (delimit_pos = strs[1].find_first_of(" \t", start_pos))
				!= string::npos) {
				if (delimit_pos > start_pos) {
					string h_id =
						strs[1].substr(start_pos, delimit_pos-start_pos-1);
					bool h_sense =
						(strs[1][delimit_pos-1] == '+'? false : true);

					g_edges_irregular.insert(pair<ContigNode, ContigNode>(
						ContigNode(t_id, false),
						ContigNode(h_id, h_sense) ));

					g_edges_irregular.insert(pair<ContigNode, ContigNode>(
						ContigNode(h_id, !h_sense),
						ContigNode(t_id, true) ));
				}
				start_pos = delimit_pos + 1;
			}
		}
		if (!strs[2].empty()) {
			start_pos = 0;
			while (start_pos < strs[2].length()
				&& (delimit_pos = strs[2].find_first_of(' ', start_pos))
				!= string::npos) {
				if (delimit_pos > start_pos) {
					string h_id =
						strs[2].substr(start_pos, delimit_pos-start_pos-1);
					bool h_sense =
						(strs[2][delimit_pos-1] == '+'? false : true);

					g_edges_irregular.insert(pair<ContigNode, ContigNode>(
						ContigNode(h_id, h_sense),
						ContigNode(t_id, false) ));

					g_edges_irregular.insert(pair<ContigNode, ContigNode>(
						ContigNode(t_id, true),
						ContigNode(h_id, !h_sense) ));

				}
				start_pos = delimit_pos + 1;
			}
		}
	}
	in_adj.close();
	if (opt::verbose > 1) {
		cerr << "g_edges_irregular:\n";
		set< pair<ContigNode, ContigNode> >::const_iterator it =
			g_edges_irregular.begin();
		for (; it != g_edges_irregular.end(); it++)
			cerr << it->first << " to " << it->second << '\n';
	}
}
#endif

/** Output a new contig. */
static ContigID outputNewContig(
	const vector<Path>& solutions,
	size_t longestPrefix, size_t longestSuffix,
	const Sequence& seq, const unsigned coverage,
	ofstream& out)
{
	ContigID id(ContigID::create());
	out << '>' << id.str() << ' ' << seq.length()
		<< ' ' << coverage << ' ';
	for (vector<Path>::const_iterator it = solutions.begin();
			it != solutions.end(); it++) {
		if (it != solutions.begin())
			out << ';';
		const ContigPath& path = *it;
		ContigPath::const_iterator
			first = path.begin() + longestPrefix,
			last = path.end() - longestSuffix;
		assert(first <= last);
		if (first < last) {
			copy(first, last - 1,
					ostream_iterator<ContigNode>(out, ","));
			out << *(last - 1);
		} else
			out << '*';
	}
	out << '\n' << seq << '\n';
	return id;
}

/** Return a consensus sequence of a and b.
 * @return an empty string if a consensus could not be found
 */
static string createConsensus(const Sequence& a, const Sequence& b)
{
	assert(a.length() == b.length());
	if (a == b)
		return a;
	string s;
	s.reserve(a.length());
	for (string::const_iterator ita = a.begin(), itb = b.begin();
			ita != a.end(); ++ita, ++itb) {
		bool mask = islower(*ita) || islower(*itb);
		char ca = toupper(*ita), cb = toupper(*itb);
		char c = ca == cb ? ca
			: ca == 'N' ? cb
			: cb == 'N' ? ca
			: 'x';
		if (c == 'x')
			return string("");
		s += mask ? tolower(c) : c;
	}
	return s;
}

/** Merge the specified two contigs,
 * With modified Overlap, we no longer assume overlap is k-1, instead,
 * do a smith-waterman alignment between contigs to determine their
 * overlap.
*/
static overlap_align mergeContigs_SW(Sequence& seq, const Sequence& s,
	const ContigNode& node, const Path& path)
{
	vector<overlap_align> overlaps;
	overlaps.reserve(1); /* only need the longest overlap */
	do {
		unsigned len = min(seq.length(), s.length());
		assert(len > opt::k - 1);
		GetLocalAlignment(seq.substr(seq.length()-len),
			s.substr(0, len),
			seq.length()-len, overlaps, false,
			opt::verbose > 1);
	} while (overlaps.empty() && chomp(seq, 'n'));

	if (overlaps.empty()) {
		cerr << "warning: the head of `" << node << "' "
			"does not match the tail of the previous contig\n"
			<< seq << '\n' << s << '\n' << path << endl;
		seq += 'n';
		seq += s;
		return overlap_align();
	} else {
		seq.erase(overlaps[0].overlap_t_pos);
		seq += overlaps[0].overlap_str;
		seq += Sequence(s, overlaps[0].overlap_h_pos+1);
		return overlaps[0];
	}
}

/** Merge the specified two contigs, default overlap is k-1,
 * generate a consensus sequence of the overlapping region. The result
 * is stored in the first argument.
 */
static void mergeContigs(Sequence& seq, const Sequence& s,
	const ContigNode& node, const Path& path)
{
	unsigned overlap = opt::k - 1;
	assert(s.length() > overlap);
	Sequence ao;
	Sequence bo(s, 0, overlap);
	Sequence o;
	do {
		assert(seq.length() > overlap);
		ao = seq.substr(seq.length() - overlap);
		o = createConsensus(ao, bo);
	} while (o.empty() && chomp(seq, 'n'));
	if (o.empty()) {
		cerr << "warning: the head of `" << node << "' "
			"does not match the tail of the previous contig\n"
			<< ao << '\n' << bo << '\n' << path << endl;
		seq += 'n';
		seq += s;
	} else {
		seq.resize(seq.length() - overlap);
		seq += o;
		seq += Sequence(s, overlap);
	}
}

static Sequence mergePath(const Path& path)
{
	Sequence seq;
	Path::const_iterator prev_it;
	for (Path::const_iterator it = path.begin();
			it != path.end(); ++it) {
		if (seq.empty()) {
			seq = getSequence(*it);
		} else {
			if (g_edges_irregular.find(
					pair<ContigNode, ContigNode>(
					*prev_it, *it))
					== g_edges_irregular.end())
				mergeContigs(seq, getSequence(*it),
					*it, path);
			else
				mergeContigs_SW(seq, getSequence(*it),
					*it, path);
		}
		prev_it = it;
	}
	return seq;
}

#if 0
/** Validate path coverage, at 95% confidence interval. */
static bool validCoverage(unsigned pathLen, unsigned pathCover)
{
	double cover_mean
		= (double)pathCover / (pathLen + opt::k - 1);
	double cover_deviation
		= sqrt(g_coverage_variance / (pathLen + opt::k - 1));
	return cover_mean <= g_coverage_mean + 1.96*cover_deviation
		&& cover_mean >= g_coverage_mean - 1.96*cover_deviation;
}
#endif

/** Calculate the ContigProperties of a path. */
static ContigProperties calculatePathProperties(const Graph& g,
		const ContigPath& path)
{
	assert(!path.empty());
	ContigPath::const_iterator first = path.begin();
	return std::accumulate(first + 1, path.end(),
			get(vertex_bundle, g, *first), AddVertexProp<Graph>(g));
}

/* Resolve ambiguous region using pairwise alignment (needleman-wunsch)
 * ('solutions' contain exactly two paths, from a source contig
 * to a dest contig)
 */
static ContigPath alignPair(const Graph& g,
		const ContigPaths& solutions, ofstream& out)
{
	assert(solutions.size() == 2);
	assert(solutions[0].size() > 1);
	assert(solutions[1].size() > 1);
	assert(solutions[0].front() == solutions[1].front());
	assert(solutions[0].back() == solutions[1].back());
	ContigPath fstSol(solutions[0].begin()+1, solutions[0].end()-1);
	ContigPath sndSol(solutions[1].begin()+1, solutions[1].end()-1);

	if (fstSol.empty() || sndSol.empty()) {
		// This entire sequence may be deleted.
		const ContigPath& sol(fstSol.empty() ? sndSol : fstSol);
		assert(!sol.empty());
		Sequence consensus(mergePath(sol));
		assert(consensus.size() > opt::k - 1);
		string::iterator first = consensus.begin() + opt::k - 1;
		transform(first, consensus.end(), first, ::tolower);
		unsigned coverage = calculatePathProperties(g, sol).coverage;
		ContigID id
			= outputNewContig(solutions, 1, 1, consensus, coverage, out);
		ContigPath path;
		path.push_back(solutions.front().front());
		path.push_back(ContigNode(id, false));
		path.push_back(solutions.front().back());
		return path;
	}

	Sequence fstPathContig(mergePath(fstSol));
	Sequence sndPathContig(mergePath(sndSol));
	if (fstPathContig == sndPathContig) {
		// A perfect match must be caused by palindrome.
		assert(fstSol.size() == 1);
		assert(sndSol.size() == 1);
		assert(fstSol.front() == ~sndSol.front());
		if (opt::verbose > 0)
			cerr << "Palindrome: " << fstSol.front().id() << '\n';
		return solutions[0];
	}

	NWAlignment align;
	unsigned match = GetGlobalAlignment(fstPathContig, sndPathContig,
		   	align);
	float identity = (float)match / align.size();
	if (opt::verbose > 1)
		cerr << align;
	if (opt::verbose > 0)
		cerr << identity
			<< (identity < opt::identity ? " (too low)\n" : "\n");
	if (identity < opt::identity)
		return ContigPath();

	unsigned coverage = calculatePathProperties(g, fstSol).coverage
		+ calculatePathProperties(g, sndSol).coverage;
	ContigID id = outputNewContig(solutions, 1, 1,
			align.consensus(), coverage, out);
	ContigPath path;
	path.push_back(solutions.front().front());
	path.push_back(ContigNode(id, false));
	path.push_back(solutions.front().back());
	return path;
}

/* Resolve ambiguous region using multiple alignment of all paths in
 * `solutions'.
 */
static ContigPath alignMulti(const Graph& g,
		const vector<Path>& solutions, ofstream& out)
{
	// Find the size of the smallest path.
	const Path& firstSol = solutions.front();
	size_t min_len = firstSol.size();
	for (vector<Path>::const_iterator it = solutions.begin() + 1;
			it != solutions.end(); ++it)
		min_len = min(min_len, it->size());

	// Find the longest prefix.
	Path vppath;
	size_t longestPrefix;
	bool commonPrefix = true;
	for (longestPrefix = 0; longestPrefix < min_len; longestPrefix++) {
		const ContigNode& common_path_node = firstSol[longestPrefix];
		for (vector<Path>::const_iterator solIter = solutions.begin();
				solIter != solutions.end(); ++solIter) {
			const ContigNode& pathnode = (*solIter)[longestPrefix];
			if (pathnode != common_path_node) {
				// Found the longest prefix.
				commonPrefix = false;
				break;
			}
		}
		if (!commonPrefix)
			break;
		vppath.push_back(common_path_node);
	}

	// Find the longest suffix.
	Path vspath;
	size_t longestSuffix;
	bool commonSuffix = true;
	for (longestSuffix = 0;	longestSuffix < min_len-longestPrefix;
			longestSuffix++) {
		const ContigNode& common_path_node
			= firstSol[firstSol.size()-longestSuffix-1];
		for (vector<Path>::const_iterator solIter = solutions.begin();
				solIter != solutions.end(); ++solIter) {
			const ContigNode& pathnode
				= (*solIter)[solIter->size()-longestSuffix-1];
			if (pathnode != common_path_node) {
				// Found the longest suffix.
				commonSuffix = false;
				break;
			}
		}
		if (!commonSuffix)
			break;
		vspath.push_back(common_path_node);
	}
	reverse(vspath.begin(), vspath.end());

	if (opt::verbose > 0 && vppath.size() + vspath.size() > 2)
		cerr << vppath << " * " << vspath << '\n';

	// Get sequence of ambiguous region in paths
	assert(longestPrefix > 0 && longestSuffix > 0);
	vector<Sequence> amb_seqs;
	unsigned coverage = 0;
	for (vector<Path>::const_iterator solIter = solutions.begin();
			solIter != solutions.end(); solIter++) {
		assert(longestPrefix + longestSuffix <= solIter->size());
		Path path(solIter->begin() + longestPrefix,
				solIter->end() - longestSuffix);
		if (!path.empty()) {
			amb_seqs.push_back(mergePath(path));
			coverage += calculatePathProperties(g, path).coverage;
		} else {
			// The prefix and suffix paths overlap by k-1 bp.
			Sequence s = getSequence(solutions[0][longestPrefix-1]);
			amb_seqs.push_back(s.substr(s.length() - opt::k + 1));
		}
	}

	unsigned matches;
	Sequence consensus = dialign(amb_seqs, matches);
	if (opt::verbose > 1)
	   	cerr << consensus << '\n';
	float identity = (float)matches / consensus.size();
	if (opt::verbose > 0)
		cerr << identity
			<< (identity < opt::identity ? " (too low)\n" : "\n");
	if (identity < opt::identity)
		return ContigPath();

	ContigID id = outputNewContig(solutions,
		longestPrefix, longestSuffix, consensus, coverage, out);
	ContigPath path(vppath);
	path.push_back(ContigNode(id, false));
	path.insert(path.end(), vspath.begin(), vspath.end());
	return path;
}

/** Align the sequences of the specified paths.
 * @return the consensus sequence
 */
static ContigPath align(const Graph& g, const vector<Path>& sequences,
		ofstream& out)
{
	assert(!sequences.empty());
	switch (sequences.size()) {
	  case 1:
		return sequences.front();
	  case 2:
		return alignPair(g, sequences, out);
	  default:
		return alignMulti(g, sequences, out);
	}
}

#if 0
static void calculateCoverageStatistics(const Graph& g)
{
	double coverage = 0, variance = 0;
	unsigned num = g.num_vertices() / 2;
	for (unsigned i = 0; i < num; ++i) {
		ContigNode u(i, false);
		double c = (double)g[u].coverage / (g[u].length - opt::k + 1);
		coverage += c;
		variance += c * c;
	}
	g_coverage_mean = coverage / num;
	g_coverage_variance = (variance - 2 * g_coverage_mean * coverage)
		/ num + g_coverage_mean * g_coverage_mean;
}
#endif

/**
 * main program routine
 */
int main(int argc, char **argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
			shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case '?': die = true; break;
		case 'k': arg >> opt::k; break;
		case 'o': arg >> opt::out; break;
		case 'p': arg >> opt::identity; break;
		case 'a': arg >> opt::numBranches; break;
		case 's': arg >> opt::consensusPath; break;
		case 'D': arg >> opt::dialign_debug; break;
		case 'M': arg >> opt::dialign_score; break;
		case 'P': arg >> opt::dialign_prob; break;
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

	if (opt::out.empty()) {
		cerr << PROGRAM ": " << "missing -o,--out option\n";
		die = true;
	}

	if (opt::consensusPath.empty()) {
		cerr << PROGRAM ": " << "missing -s,--consensus option\n";
		die = true;
	}

	if (argc - optind < 3) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	const char *contigFile = argv[optind++];
	string adjFile(argv[optind++]);
	string allPaths(argv[optind++]);

	// Load the graph from the adjacency file
	Graph g;
	ifstream fin(adjFile.c_str());
	fin >> g;
	assert(fin.eof());

	// Read contigs
	Contigs& contigs = g_contigs;
	{
		FastaReader in(contigFile, FastaReader::NO_FOLD_CASE);
		for (FastaRecord rec; in >> rec;) {
			ContigID id(rec.id);
			assert(contigs.size() == id);
			contigs.push_back(rec.seq);
		}
		assert(in.eof());
		assert(!contigs.empty());
		opt::colourSpace = isdigit(contigs[0][0]);
	}
	ContigID::lock();

#if 0
	// Get contig k-mer-coverage statistics.
	calculateCoverageStatistics(g);
#endif

	vector<string> pathIDs;
	vector<bool> isAmbPath;
	ContigPaths paths = readPath(allPaths, pathIDs, isAmbPath);
	stats.numAmbPaths = g_ambpath_contig.size();
	if (opt::verbose > 0)
		cerr << "Read " << paths.size() << " paths\n";

	// Start numbering new contigs from the last
	ContigID::setNextContigID(pathIDs.back());

	// Prepare output fasta file
	ofstream fa(opt::consensusPath.c_str());
	assert(fa.good());

	init_parameters();
	set_parameters_dna();
	para->DEBUG = opt::dialign_debug;
	para->SCR_MATRIX_FILE_NAME = (char*)opt::dialign_score.c_str();
	para->DIAG_PROB_FILE_NAME = (char*)opt::dialign_prob.c_str();
	initDialign();

	// Contigs that were seen in a consensus.
	vector<bool> seen(contigs.size());

	// resolve ambiguous paths recorded in g_ambpath_contig
	for (AmbPath2Contig::iterator ambIt = g_ambpath_contig.begin();
			ambIt != g_ambpath_contig.end(); ambIt++) {
		AmbPathConstraint apConstraint = ambIt->first;
		if (opt::verbose > 0)
			cerr << "\n* " << apConstraint.source << ' '
				<< apConstraint.dist << "N "
				<< apConstraint.dest << '\n';

		Constraints constraints;
		constraints.push_back(Constraint(apConstraint.dest,
			apConstraint.dist + opt::allowedError));

		ContigPaths solutions;
		unsigned numVisited = 0;
		constrainedSearch(g, apConstraint.source,
			constraints, solutions, numVisited);
		bool tooComplex = numVisited >= opt::maxCost;

		for (ContigPaths::iterator solIt = solutions.begin();
				solIt != solutions.end(); solIt++)
			solIt->insert(solIt->begin(), apConstraint.source);

		bool tooManySolutions = solutions.size() > opt::numBranches;
		if (tooComplex) {
			stats.tooComplex++;
			if (opt::verbose > 0)
				cerr << solutions.size() << " paths (too complex)\n";
		} else if (tooManySolutions) {
			stats.numTooManySolutions++;
			if (opt::verbose > 0)
				cerr << solutions.size() << " paths (too many)\n";
		} else if (solutions.empty()) {
			stats.numNoSolutions++;
			if (opt::verbose > 0)
				cerr << "no paths\n";
		} else {
			if (opt::verbose > 1)
				copy(solutions.begin(), solutions.end(),
						ostream_iterator<ContigPath>(cerr, "\n"));
			ambIt->second = align(g, solutions, fa);
			if (!ambIt->second.empty()) {
				stats.numMerged++;
				if (solutions.size() > 1) {
					// Mark contigs that are used in a consensus.
					markSeen(seen, solutions, true);
				}
				if (opt::verbose > 0)
					cerr << ambIt->second << '\n';
			} else
				stats.notMerged++;
		}
	}
	if (opt::verbose > 0)
		cerr << '\n';

	// Unmark contigs that are used in a path.
	for (AmbPath2Contig::iterator it = g_ambpath_contig.begin();
			it != g_ambpath_contig.end(); it++)
		markSeen(seen, it->second, false);
	markSeen(seen, paths, false);

	ofstream out(opt::out.c_str());
	assert(out.good());

	// Output those contigs that were not seen in ambiguous path.
	for (unsigned id = 0; id < contigs.size(); ++id)
		if (seen[id])
			out << ContigID(id) << '\n';

	for (ContigPaths::const_iterator path = paths.begin();
			path != paths.end(); ++path) {
		unsigned i = path - paths.begin();
		if (!isAmbPath[i]) {
			out << pathIDs[i] << '\t' << *path << '\n';
			continue;
		}

		assert(path->size() > 2);
		Path cur_path;
		cur_path.push_back(path->front());
		for (Path::const_iterator prev = path->begin(),
				it = path->begin() + 1, next = path->begin() + 2;
				it != path->end(); ++prev, ++it, ++next) {
			if (!it->ambiguous()) {
				cur_path.push_back(*it);
				continue;
			}

			//replace Ns by resolved new contig
			assert(next != path->end());
			AmbPath2Contig::iterator ambIt = g_ambpath_contig.find(
				AmbPathConstraint(*prev, *next, -it->id()));
			assert(ambIt != g_ambpath_contig.end());
			const ContigPath& solution = ambIt->second;
			if (!solution.empty()) {
				assert(solution.size() > 1);
				cur_path.insert(cur_path.end(),
						solution.begin() + 1, solution.end() - 1);
			} else
				cur_path.push_back(*it);
		}
		out << pathIDs[i] << '\t' << cur_path << '\n';
	}

	free_prob_dist(pdist);
	free(para);

	out.close();
	fa.close();

	cerr <<
		"Ambiguous paths: " << stats.numAmbPaths << "\n"
		"Merged:          " << stats.numMerged << "\n"
		"No paths:        " << stats.numNoSolutions << "\n"
		"Too many paths:  " << stats.numTooManySolutions << "\n"
		"Too complex:     " << stats.tooComplex << "\n"
		"Dissimilar:      " << stats.notMerged << "\n";
}
