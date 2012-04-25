#include "config.h"
#include "Common/Options.h"
#include "ContigNode.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "DataLayer/Options.h"
#include "Dictionary.h"
#include "FastaReader.h"
#include "Histogram.h"
#include "IOUtil.h"
#include "MemoryUtil.h"
#include "smith_waterman.h"
#include "Sequence.h"
#include "StringUtil.h"
#include "Uncompress.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include "Graph/GraphUtil.h"
#include "Graph/Options.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

#define PROGRAM "MergeContigs"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Shaun Jackman.\n"
"\n"
"Copyright 2012 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... FASTA [OVERLAP] PATH\n"
"Merge paths of contigs to create larger contigs.\n"
"  FASTA    contigs in FASTA format\n"
"  OVERLAP  contig overlap graph\n"
"  PATH     sequences of contig IDs\n"
"\n"
"  -k, --kmer=KMER_SIZE  k-mer size\n"
"  -o, --out=FILE        output the merged contigs to FILE [stdout]\n"
"  -g, --graph=FILE      write the contig overlap graph to FILE\n"
"      --merged          output only merged contigs\n"
"      --adj             output the graph in adj format\n"
"      --dot             output the graph in dot format [default]\n"
"      --dot-meancov     same as above but give the mean coverage\n"
"      --sam             output the graph in SAM format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k; // used by ContigProperties

	/** Output FASTA path. */
	static string out = "-";

	/** Output graph path. */
	static string graphPath;

	/** Output graph format. */
	int format = DOT;

	/** Output only merged contigs. */
	int onlyMerged;

	/** Minimum overlap. */
	static unsigned minOverlap = 20;

	/** Minimum alignment identity. */
	static float minIdentity = 0.9;
}

static const char shortopts[] = "g:k:o:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "adj", no_argument, &opt::format, ADJ },
	{ "dot", no_argument, &opt::format, DOT },
	{ "dot-meancov", no_argument, &opt::format, DOT_MEANCOV },
	{ "sam", no_argument, &opt::format, SAM },
	{ "graph", required_argument, NULL, 'g' },
	{ "kmer",        required_argument, NULL, 'k' },
	{ "merged",      no_argument,       &opt::onlyMerged, 1 },
	{ "out",         required_argument, NULL, 'o' },
	{ "path",        required_argument, NULL, 'p' },
	{ "verbose",     no_argument,       NULL, 'v' },
	{ "help",        no_argument,       NULL, OPT_HELP },
	{ "version",     no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/* A contig sequence. */
struct Contig {
	Contig(const string& comment, const string& seq)
		: comment(comment), seq(seq) { }
	Contig(const FastaRecord& o) : comment(o.comment), seq(o.seq) { }
	string comment;
	string seq;
};

/** The contig sequences. */
typedef vector<Contig> Contigs;

/** Return the sequence of the specified contig node. The sequence
 * may be ambiguous or reverse complemented.
 */
static Sequence sequence(const Contigs& contigs, const ContigNode& id)
{
	if (id.ambiguous()) {
		string s(id.ambiguousSequence());
		if (s.length() < opt::k)
			transform(s.begin(), s.end(), s.begin(), ::tolower);
		return string(opt::k - 1, 'N') + s;
	} else {
		const Sequence& seq = contigs[id.id()].seq;
		return id.sense() ? reverseComplement(seq) : seq;
	}
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
			: ambiguityIsSubset(ca, cb) ? ambiguityOr(ca, cb)
			: 'x';
		if (c == 'x')
			return string("");
		s += mask ? tolower(c) : c;
	}
	return s;
}

typedef ContigGraph<DirectedGraph<ContigProperties, Distance> > Graph;
typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;

/** Return the properties of the specified vertex, unless u is
 * ambiguous, in which case return the length of the ambiguous
 * sequence.
 */
static inline
ContigProperties get(vertex_bundle_t, const Graph& g, ContigNode u)
{
	return u.ambiguous()
		? ContigProperties(u.length() + opt::k - 1, 0)
		: g[u];
}

/** Append the sequence of contig v to seq. */
static void mergeContigs(const Graph& g, const Contigs& contigs,
		vertex_descriptor u, vertex_descriptor v,
		Sequence& seq, const ContigPath& path)
{
	int d = get(edge_bundle, g, u, v).distance;
	assert(d < 0);
	unsigned overlap = -d;
	const Sequence& s = sequence(contigs, v);
	assert(s.length() > overlap);
	Sequence ao;
	Sequence bo(s, 0, overlap);
	Sequence o;
	do {
		assert(seq.length() > overlap);
		ao = seq.substr(seq.length() - overlap);
		o = createConsensus(ao, bo);
		if (!o.empty()) {
			seq.resize(seq.length() - overlap);
			seq += o;
			seq += Sequence(s, overlap);
			return;
		}
	} while (chomp(seq, 'n'));

	// Try an overlap alignment.
	if (opt::verbose > 2)
		cerr << '\n';
	vector<overlap_align> overlaps;
	alignOverlap(ao, bo, 0, overlaps, false, opt::verbose > 2);
	bool good = false;
	if (!overlaps.empty()) {
		assert(overlaps.size() == 1);
		const overlap_align& o = overlaps.front();
		unsigned matches = o.overlap_match;
		const string& consensus = o.overlap_str;
		float identity = (float)matches / consensus.size();
		good = matches >= opt::minOverlap
			&& identity >= opt::minIdentity;
		if (opt::verbose > 2)
			cerr << matches << " / " << consensus.size()
				<< " = " << identity
				<< (matches < opt::minOverlap ? " (too few)"
						: identity < opt::minIdentity ? " (too low)"
						: " (good)") << '\n';
	}
	if (good) {
		assert(overlaps.size() == 1);
		const overlap_align& o = overlaps.front();
		seq.erase(seq.length() - overlap + o.overlap_t_pos);
		seq += o.overlap_str;
		seq += Sequence(s, o.overlap_h_pos + 1);
	} else {
		cerr << "warning: the head of " << get(vertex_name, g, v)
			<< " does not match the tail of the previous contig\n"
			<< ao << '\n' << bo << '\n' << path << endl;
		seq += 'n';
		seq += s;
	}
}

/** Return a FASTA comment for the specified path. */
static void pathToComment(ostream& out,
		const Graph& g, const ContigPath& path)
{
	assert(path.size() > 1);
	out << get(vertex_name, g, path.front());
	if (path.size() == 3)
		out << ',' << get(vertex_name, g, path[1]);
	else if (path.size() > 3)
		out << ",...";
	out << ',' << get(vertex_name, g, path.back());
}

/** Merge the specified path. */
static Contig mergePath(const Graph& g, const Contigs& contigs,
		const ContigPath& path)
{
	Sequence seq;
	unsigned coverage = 0;
	for (ContigPath::const_iterator it = path.begin();
			it != path.end(); ++it) {
		if (!it->ambiguous())
			coverage += g[*it].coverage;
		if (seq.empty()) {
			seq = sequence(contigs, *it);
		} else {
			assert(it != path.begin());
			mergeContigs(g, contigs, *(it-1), *it, seq, path);
		}
	}
	ostringstream ss;
	ss << seq.size() << ' ' << coverage << ' ';
	pathToComment(ss, g, path);
	return Contig(ss.str(), seq);
}

/** A container of ContigPath. */
typedef vector<ContigPath> ContigPaths;

/** Read contig paths from the specified file.
 * @param ids [out] the string ID of the paths
 */
static ContigPaths readPaths(const string& inPath,
		vector<string>* ids = NULL)
{
	if (ids != NULL)
		assert(ids->empty());
	ifstream fin(inPath.c_str());
	if (opt::verbose > 0)
		cerr << "Reading `" << inPath << "'..." << endl;
	if (inPath != "-")
		assert_good(fin, inPath);
	istream& in = inPath == "-" ? cin : fin;

	unsigned count = 0;
	ContigPaths paths;
	string id;
	ContigPath path;
	while (in >> id >> path) {
		paths.push_back(path);
		if (ids != NULL)
			ids->push_back(id);

		++count;
		if (opt::verbose > 1 && count % 1000000 == 0)
			cerr << "Read " << count << " paths. "
				"Using " << toSI(getMemoryUsage())
				<< "B of memory.\n";
	}
	if (opt::verbose > 0)
		cerr << "Read " << count << " paths. "
			"Using " << toSI(getMemoryUsage()) << "B of memory.\n";
	assert(in.eof());
	return paths;
}

/** Finds all contigs used in each path in paths, and
 * marks them as seen in the vector seen. */
static void seenContigs(vector<bool>& seen, const ContigPaths& paths)
{
	for (ContigPaths::const_iterator it = paths.begin();
			it != paths.end(); ++it)
		for (ContigPath::const_iterator itc = it->begin();
				itc != it->end(); ++itc)
			if (itc->id() < seen.size())
				seen[itc->id()] = true;
}

/** Mark contigs for removal. An empty path indicates that a contig
 * should be removed.
 */
static void markRemovedContigs(vector<bool>& marked,
		const vector<string>& pathIDs, const ContigPaths& paths)
{
	for (ContigPaths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		if (it->empty()) {
			size_t i = get(g_contigNames,
					pathIDs[it - paths.begin()]);
			assert(i < marked.size());
			marked[i] = true;
		}
	}
}

/** Output the updated overlap graph. */
static void outputGraph(Graph& g,
		const vector<string>& pathIDs, const ContigPaths& paths,
		const string& commandLine)
{
	typedef graph_traits<Graph>::vertex_descriptor V;

	// Add the path vertices.
	g_contigNames.unlock();
	for (ContigPaths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		const ContigPath& path = *it;
		const string& id = pathIDs[it - paths.begin()];
		if (!path.empty()) {
			V u = merge(g, path.begin(), path.end());
			put(vertex_name, g, u, id);
		}
	}
	g_contigNames.lock();

	// Remove the vertices that are used in paths.
	for (ContigPaths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		const ContigPath& path = *it;
		const string& id = pathIDs[it - paths.begin()];
		if (path.empty()) {
			remove_vertex(ContigNode(id, false), g);
		} else {
			remove_vertex_if(g, path.begin(), path.end(),
					not1(std::mem_fun_ref(&ContigNode::ambiguous)));
		}
	}

	// Output the graph.
	const string& graphPath = opt::graphPath;
	assert(!graphPath.empty());
	if (opt::verbose > 0)
		cerr << "Writing `" << graphPath << "'..." << endl;
	ofstream fout(graphPath.c_str());
	assert_good(fout, graphPath);
	write_graph(fout, g, PROGRAM, commandLine);
	assert_good(fout, graphPath);
	if (opt::verbose > 0)
		printGraphStats(cerr, g);
}

int main(int argc, char** argv)
{
	opt::trimMasked = false;

	string commandLine;
	{
		ostringstream ss;
		char** last = argv + argc - 1;
		copy(argv, last, ostream_iterator<const char *>(ss, " "));
		ss << *last;
		commandLine = ss.str();
	}

	bool die = false;
	for (int c; (c = getopt_long(argc, argv,
					shortopts, longopts, NULL)) != -1;) {
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
			case '?': die = true; break;
			case 'g': arg >> opt::graphPath; break;
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
		if (optarg != NULL && !arg.eof()) {
			cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
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

	if (argc - optind < 2) {
		cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (argc - optind > 3) {
		cerr << PROGRAM ": too many arguments\n";
		die = true;
	}

	if (die) {
		cerr << "Try `" << PROGRAM
			<< " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	const char* contigFile = argv[optind++];
	string adjPath, mergedPathFile;
	Graph g;
	if (argc - optind > 1) {
		adjPath = string(argv[optind++]);

		// Read the contig adjacency graph.
		if (opt::verbose > 0)
			cerr << "Reading `" << adjPath << "'..." << endl;
		ifstream fin(adjPath.c_str());
		assert_good(fin, adjPath);
		fin >> g;
		assert(fin.eof());
		if (opt::verbose > 0)
			cerr << "Read " << num_vertices(g) << " vertices. "
				"Using " << toSI(getMemoryUsage())
				<< "B of memory.\n";
	}
	mergedPathFile = string(argv[optind++]);

	// Read the contig sequence.
	Contigs contigs;
	{
		if (opt::verbose > 0)
			cerr << "Reading `" << contigFile << "'..." << endl;
		unsigned count = 0;
		FastaReader in(contigFile, FastaReader::NO_FOLD_CASE);
		for (FastaRecord rec; in >> rec;) {
			if (!adjPath.empty()
					&& g_contigNames.count(rec.id) == 0)
				continue;
			if (adjPath.empty())
				put(g_contigNames, contigs.size(), rec.id);
			else
				assert(get(g_contigNames, rec.id) == contigs.size());
			contigs.push_back(rec);

			++count;
			if (opt::verbose > 1 && count % 1000000 == 0)
				cerr << "Read " << count << " sequences. "
					"Using " << toSI(getMemoryUsage())
					<< "B of memory.\n";
		}
		if (opt::verbose > 0)
			cerr << "Read " << count << " sequences. "
				"Using " << toSI(getMemoryUsage())
				<< "B of memory.\n";
		assert(in.eof());
		assert(!contigs.empty());
		opt::colourSpace = isdigit(contigs[0].seq[0]);
		g_contigNames.lock();
	}

	vector<string> pathIDs;
	ContigPaths paths = readPaths(mergedPathFile, &pathIDs);

	// Record all the contigs that are in a path.
	vector<bool> seen(contigs.size());
	seenContigs(seen, paths);
	markRemovedContigs(seen, pathIDs, paths);

	// Output those contigs that were not seen in a path.
	Histogram lengthHistogram;
	ofstream fout;
	ostream& out = opt::out == "-" ? cout
		: (fout.open(opt::out.c_str()), fout);
	assert_good(out, opt::out);
	if (!opt::onlyMerged) {
		for (Contigs::const_iterator it = contigs.begin();
				it != contigs.end(); ++it) {
			ContigID id(it - contigs.begin());
			if (!seen[id]) {
				const Contig& contig = *it;
				out << '>' << get(g_contigNames, id);
				if (!contig.comment.empty())
					out << ' ' << contig.comment;
				out << '\n' << contig.seq << '\n';
				if (opt::verbose > 0)
					lengthHistogram.insert(
						count_if(contig.seq.begin(), contig.seq.end(),
							isACGT));
			}
		}
	}

	if (adjPath.empty())
		return 0;

	unsigned npaths = 0;
	for (ContigPaths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		const ContigPath& path = *it;
		if (path.empty())
			continue;
		Contig contig = mergePath(g, contigs, path);
		out << '>' << pathIDs[it - paths.begin()]
			<< ' ' << contig.comment << '\n'
			<< contig.seq << '\n';
		assert_good(out, opt::out);
		npaths++;
		if (opt::verbose > 0)
			lengthHistogram.insert(
					count_if(contig.seq.begin(), contig.seq.end(),
						isACGT));
	}

	if (npaths == 0)
		return 0;

	float minCov = numeric_limits<float>::infinity(),
		minCovUsed = numeric_limits<float>::infinity();
	for (unsigned i = 0; i < contigs.size(); i++) {
		ContigProperties vp = g[ContigNode(i, false)];
		if (vp.coverage == 0)
			continue;
		assert((int)vp.length - opt::k + 1 > 0);
		float cov = (float)vp.coverage / (vp.length - opt::k + 1);
		minCov = min(minCov, cov);
		if (seen[i])
			minCovUsed = min(minCovUsed, cov);
	}

	if (!opt::graphPath.empty())
		outputGraph(g, pathIDs, paths, commandLine);

	cerr << "The minimum coverage of single-end contigs is "
		<< minCov << ".\n"
		<< "The minimum coverage of merged contigs is "
		<< minCovUsed << ".\n";
	if (minCov < minCovUsed)
		cerr << "Consider increasing the coverage threshold "
			"parameter, c, to " << minCovUsed << ".\n";

	if (opt::verbose > 0) {
		const unsigned STATS_MIN_LENGTH = 200; // bp
		printContiguityStats(cerr, lengthHistogram, STATS_MIN_LENGTH)
			<< '\t' << opt::out << '\n';
	}
	return 0;
}
