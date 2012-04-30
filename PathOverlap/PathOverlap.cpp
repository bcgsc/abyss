#include "config.h"
#include "ContigID.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "Functional.h"
#include "IOUtil.h"
#include "Uncompress.h"
#include "Graph/ContigGraph.h"
#include "Graph/ContigGraphAlgorithms.h"
#include "Graph/DirectedGraph.h"
#include "Graph/GraphIO.h"
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstring> // for strerror
#include <cstdlib>
#include <functional>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <getopt.h>
#include <map>
#include <vector>

using namespace std;

#define PROGRAM "PathOverlap"

static const char *VERSION_MESSAGE =
PROGRAM " (ABySS) " VERSION "\n"
"Written by Shaun Jackman and Tony Raymond.\n"
"\n"
"Copyright 2012 Canada's Michael Smith Genome Science Centre\n";

static const char *USAGE_MESSAGE =
"Usage: " PROGRAM " [OPTION]... ADJ PATH\n"
"Find paths that overlap. Either output the graph of overlapping\n"
"paths, assemble overlapping paths into larger paths, or trim the\n"
"overlapping paths.\n"
"  ADJ   contig adjacency graph\n"
"  PATH  sequences of contig IDs\n"
"\n"
"  -k, --kmer=N          k-mer size\n"
"  -g, --graph=FILE      write the contig adjacency graph to FILE\n"
"  -r, --repeats=FILE    write repeat contigs to FILE\n"
"      --overlap         find overlapping paths [default]\n"
"      --assemble        assemble overlapping paths\n"
"      --trim            trim overlapping paths\n"
"      --adj             output the graph in adj format [default]\n"
"      --dot             output the graph in dot format\n"
"      --sam             output the graph in SAM format\n"
"  -v, --verbose         display verbose output\n"
"      --help            display this help and exit\n"
"      --version         output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {
	unsigned k;

	/** Output format. */
	int format; // used by ContigProperties

	/** Write the contig adjacency graph to this file. */
	static string graphPath;

	/** Output the IDs of contigs in overlaps to this file. */
	static string repeatContigs;

	/** Mode of operation. */
	enum {
		/** Find overlapping paths, do not assemble. */
		OVERLAP,
		/** Assemble overlapping paths. */
		ASSEMBLE,
		/** Trim overlapping paths. */
		TRIM,
	};
	static int mode;

	static int verbose;
}

static const char* shortopts = "g:k:r:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "graph",        required_argument, NULL, 'g' },
	{ "kmer",         required_argument, NULL, 'k' },
	{ "assemble",     no_argument,       &opt::mode, opt::ASSEMBLE },
	{ "overlap",      no_argument,       &opt::mode, opt::OVERLAP },
	{ "trim",         no_argument,       &opt::mode, opt::TRIM },
	{ "adj",          no_argument,       &opt::format, ADJ, },
	{ "dot",          no_argument,       &opt::format, DOT, },
	{ "sam",          no_argument,       &opt::format, SAM, },
	{ "repeats",      required_argument, NULL, 'r' },
	{ "verbose",      no_argument,       NULL, 'v' },
	{ "help",         no_argument,       NULL, OPT_HELP },
	{ "version",      no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

/** A vertex of the overlap graph. */
struct Vertex {
	unsigned id;
	bool sense;

	/** The number of single-end contigs. */
	static unsigned s_offset;

	Vertex(unsigned id, bool sense)
		: id(id), sense(sense) { }

	bool operator ==(const Vertex& v) const
	{
		return id == v.id && sense == v.sense;
	}

	operator ContigNode() const
	{
		return ContigNode(s_offset + id, sense);
	}
};

unsigned Vertex::s_offset;

/** An alignment of two overlapping contigs. */
struct Overlap {
	Vertex source;
	Vertex target;

	/** Overlap measured in number of contigs. */
	unsigned overlap;

	/** Overlap measured in bp. */
	int distance;

	Overlap(const Vertex& source, const Vertex& target,
			unsigned overlap, int distance)
		: source(source), target(target),
		overlap(overlap), distance(distance) { }
};

/** The contig IDs that have been removed from paths. */
static vector<ContigID> s_trimmedContigs;

/** The contig graph. */
typedef DirectedGraph<ContigProperties, Distance> DG;
typedef ContigGraph<DG> Graph;

typedef vector<ContigPath> Paths;

/** Return whether this vertex is a path or a contig. */
static bool isPath(const ContigNode& u)
{
	return u.id() >= Vertex::s_offset;
}

/** Return a path, complemented if necessary. */
static ContigPath getPath(const Paths& paths, const ContigNode& u)
{
	if (isPath(u)) {
		unsigned i = u.id() - Vertex::s_offset;
		return u.sense() ? reverseComplement(paths[i]) : paths[i];
	} else
		return ContigPath(1, u);
}

/** Read contig paths from the specified file.
 * @param g the contig adjacency graph
 * @param inPath the file of contig paths
 * @param[out] pathIDs the path IDs
 * @return the paths
 */
static Paths readPaths(Graph& g,
		const string& inPath, vector<string>& pathIDs)
{
	assert(pathIDs.empty());
	ifstream fin(inPath.c_str());
	if (opt::verbose > 0)
		cerr << "Reading `" << inPath << "'..." << endl;
	if (inPath != "-")
		assert_good(fin, inPath);
	istream& in = inPath == "-" ? cin : fin;

	assert_good(in, inPath);
	Paths paths;
	string id;
	ContigPath path;
	while (in >> id >> path) {
		if (path.empty()) {
			// Remove this contig from the graph.
			ContigNode u(id, false);
			clear_vertex(u, g);
			remove_vertex(u, g);
		} else {
			pathIDs.push_back(id);
			paths.push_back(path);
		}
	}
	assert(in.eof());
	return paths;
}

typedef multimap<ContigNode, Vertex> SeedMap;

/** Index the first and last contig of each path to facilitate finding
 * overlaps between paths. */
static SeedMap makeSeedMap(const Paths& paths)
{
	SeedMap seedMap;
	for (Paths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		if (it->empty())
			continue;
		assert(!it->front().ambiguous());
		seedMap.insert(make_pair(it->front(),
					Vertex(it - paths.begin(), false)));
		assert(!it->back().ambiguous());
		seedMap.insert(make_pair(~it->back(),
					Vertex(it - paths.begin(), true)));
	}
	return seedMap;
}

/** Check whether path starts with the sequence [first, last). */
static bool startsWith(ContigPath path, bool rc,
		ContigPath::const_iterator first,
		ContigPath::const_iterator last)
{
	if (rc)
		reverseComplement(path.begin(), path.end());
	assert(*first == path.front());
	assert(first < last);
	return unsigned(last - first) > path.size() ? false
		: equal(first, last, path.begin());
}

/** Check whether path starts with the sequence [first, last). */
static unsigned findOverlap(const Graph& g,
		const Paths& paths,
		ContigPath::const_iterator first,
		ContigPath::const_iterator last,
		const Vertex& v, int &distance)
{
	if (!startsWith(paths[v.id], v.sense, first, last))
		return 0;
	distance = -addProp(g, first, last).length;
	return last - first;
}

typedef vector<Overlap> Overlaps;

/** Find every path that overlaps with the specified path. */
static void findOverlaps(const Graph& g,
		const Paths& paths, const SeedMap& seedMap,
		const Vertex& v, Overlaps& overlaps)
{
	ContigPath rc;
	if (v.sense) {
		rc = paths[v.id];
		reverseComplement(rc.begin(), rc.end());
	}
	const ContigPath& path = v.sense ? rc : paths[v.id];

	for (ContigPath::const_iterator it = path.begin();
			it != path.end(); ++it) {
		if (it->ambiguous())
			continue;

		pair<SeedMap::const_iterator, SeedMap::const_iterator>
			range = seedMap.equal_range(*it);
		for (SeedMap::const_iterator seed = range.first;
				seed != range.second; ++seed) {
			if (v == seed->second)
				continue;
			int distance = 0;
			unsigned overlap = findOverlap(g, paths, it, path.end(),
					   seed->second, distance);
			if (overlap > 0)
				overlaps.push_back(Overlap(v, seed->second,
					overlap, distance));

		}
	}
}

/** Find every pair of overlapping paths. */
static Overlaps findOverlaps(const Graph& g, const Paths& paths)
{
	SeedMap seedMap = makeSeedMap(paths);

	Overlaps overlaps;
	for (Paths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		unsigned i = it - paths.begin();
		findOverlaps(g, paths, seedMap, Vertex(i, false), overlaps);
		findOverlaps(g, paths, seedMap, Vertex(i, true), overlaps);
	}
	return overlaps;
}

/** Record the trimmed contigs. */
static void recordTrimmedContigs(
		ContigPath::const_iterator first,
		ContigPath::const_iterator last)
{
	for (ContigPath::const_iterator it = first; it != last; ++it)
		if (!it->ambiguous())
			s_trimmedContigs.push_back(ContigID(*it));
}

/** Remove ambiguous contigs from the ends of the path. */
static void removeAmbiguousContigs(ContigPath& path)
{
	if (!path.empty() && path.back().ambiguous())
		path.erase(path.end() - 1);
	if (!path.empty() && path.front().ambiguous())
		path.erase(path.begin());
}

/** Remove the overlapping portion of the specified contig. */
static void removeContigs(ContigPath& path,
		unsigned first, unsigned last)
{
	assert(first <= path.size());
	assert(last <= path.size());
	if (first < last) {
		recordTrimmedContigs(path.begin(), path.begin() + first);
		recordTrimmedContigs(path.begin() + last, path.end());
		path.erase(path.begin() + last, path.end());
		path.erase(path.begin(), path.begin() + first);
	} else {
		recordTrimmedContigs(path.begin(), path.end());
		path.clear();
	}
	removeAmbiguousContigs(path);
}

/** Find the largest overlap for each contig and remove it. */
static void trimOverlaps(Paths& paths, const Overlaps& overlaps)
{
	vector<unsigned> removed[2];
	removed[0].resize(paths.size());
	removed[1].resize(paths.size());

	for (Overlaps::const_iterator it = overlaps.begin();
			it != overlaps.end(); ++it) {
		unsigned& a = removed[!it->source.sense][it->source.id];
		unsigned& b = removed[it->target.sense][it->target.id];
		a = max(a, it->overlap);
		b = max(b, it->overlap);
	}

	for (Paths::iterator it = paths.begin(); it != paths.end(); ++it)
		removeContigs(*it, removed[0][it - paths.begin()],
				it->size() - removed[1][it - paths.begin()]);
}

/** Trim the ends of paths that overlap another path. */
static void trimOverlaps(const Graph& g, Paths& paths)
{
	for (Overlaps overlaps = findOverlaps(g, paths);
			!overlaps.empty();
			overlaps = findOverlaps(g, paths)) {
		cerr << "Found " << overlaps.size() / 2 << " overlaps.\n";
		trimOverlaps(paths, overlaps);
	}
}

static inline
ContigProperties get(vertex_bundle_t, const Graph& g, ContigNode u)
{
	return u.ambiguous()
		? ContigProperties(u.length() + opt::k - 1, 0)
		: g[u];
}

/** Add the path overlap edges to the specified graph. */
static void addPathOverlapEdges(Graph& g,
		const Paths& paths, const vector<string>& pathIDs,
		const Overlaps& overlaps)
{
	const bool allowParallelEdge = opt::mode == opt::ASSEMBLE;

	// Add the path vertices.
	ContigID::unlock();
	for (Paths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		if (it->empty())
			continue;
		ContigID::insert(pathIDs[it - paths.begin()]);
		merge(g, it->begin(), it->end());
	}
	ContigID::lock();

	// Remove the single-end contigs that are in paths.
	for (Paths::const_iterator it = paths.begin();
			it != paths.end(); ++it)
		remove_vertex_if(g, it->begin(), it->end(),
				not1(std::mem_fun_ref(&ContigNode::ambiguous)));

	// Add the path edges.
	for (Overlaps::const_iterator it = overlaps.begin();
			it != overlaps.end(); ++it) {
		Vertex u = it->source, v = it->target;
		if (allowParallelEdge || !edge(u, v, g).second)
			add_edge(u, v, it->distance, static_cast<DG&>(g));
		else if (opt::verbose > 0)
			cerr << "ambiguous overlap: "
				<< ContigNode(u) << " -> " << ContigNode(v) << '\n';
	}
}

typedef graph_traits<Graph>::edge_descriptor edge_descriptor;

/** A property map giving the number of contigs by which two paths
 * overlap. */
typedef map<edge_descriptor, unsigned> OverlapMap;

/** Return the number of contigs by which the two paths overlap. */
static unsigned getOverlap(const OverlapMap& pmap,
		const ContigNode& u, const ContigNode& v)
{
	if (isPath(u) && isPath(v)) {
		// Both vertices are paths.
		OverlapMap::const_iterator it = pmap.find(
				edge_descriptor(u, v));
		return it == pmap.end() ? 0 : it->second;
	} else {
		// One of the two vertices is a contig.
		return 0;
	}
}

/** Merge a sequence of overlapping paths. */
static ContigPath mergePaths(const Paths& paths,
		const OverlapMap& overlaps, const ContigPath& merge)
{
	assert(!merge.empty());
	ContigNode u = merge.front();
	ContigPath path(getPath(paths, u));
	for (ContigPath::const_iterator it = merge.begin() + 1;
			it != merge.end(); ++it) {
		ContigNode v = *it;
		ContigPath vpath(getPath(paths, v));
		unsigned overlap = getOverlap(overlaps, u, v);
		assert(path.size() > overlap);
		assert(vpath.size() > overlap);
		assert(equal(path.end() - overlap, path.end(),
					vpath.begin()));
		path.insert(path.end(), vpath.begin() + overlap, vpath.end());
		u = v;
	}
	return path;
}

/** Return true if the edge e is a path overlap. */
struct IsPathOverlap : unary_function<edge_descriptor, bool> {
	IsPathOverlap(const Graph& g, const OverlapMap& pmap)
		: m_g(g), m_pmap(pmap) { }
	bool operator()(edge_descriptor e) const
	{
		return getOverlap(m_pmap, source(e, m_g), target(e, m_g));
	}
  private:
	const Graph& m_g;
	const OverlapMap& m_pmap;
};

/** Assemble overlapping paths. */
static void assembleOverlappingPaths(Graph& g,
		Paths& paths, vector<string>& pathIDs)
{
	if (paths.empty())
		return;

	// Find overlapping paths.
	Overlaps overlaps = findOverlaps(g, paths);
	addPathOverlapEdges(g, paths, pathIDs, overlaps);

	// Create a property map of path overlaps.
	OverlapMap overlapMap;
	for (Overlaps::const_iterator it = overlaps.begin();
			it != overlaps.end(); ++it)
		overlapMap.insert(OverlapMap::value_type(
				OverlapMap::key_type(it->source, it->target),
				it->overlap));

	// Assemble unambiguously overlapping paths.
	Paths merges;
	assemble_if(g, back_inserter(merges),
			IsPathOverlap(g, overlapMap));

	// Merge overlapping paths.
	assert(!pathIDs.empty());
	ContigID::setNextContigID(pathIDs.back());
	for (Paths::const_iterator it = merges.begin();
			it != merges.end(); ++it) {
		ContigID id = ContigID::create();
		if (opt::verbose > 0)
			cerr << id << '\t' << *it << '\n';
		pathIDs.push_back((string)id.str());
		paths.push_back(mergePaths(paths, overlapMap, *it));

		// Remove the merged paths.
		for (ContigPath::const_iterator it2 = it->begin();
				it2 != it->end(); ++it2) {
			if (isPath(*it2))
				paths[it2->id() - Vertex::s_offset].clear();
		}
	}
}

int main(int argc, char** argv)
{
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
			case 'r': arg >> opt::repeatContigs; break;
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

	const char *adjPath = argv[optind++];
	if (opt::verbose > 0)
		cerr << "Reading `" << adjPath << "'..." << endl;
	ifstream fin(adjPath);
	assert_good(fin, adjPath);
	Graph g;
	fin >> g;
	Vertex::s_offset = g.num_vertices() / 2;

	string pathsFile(argv[optind++]);
	vector<string> pathIDs;
	Paths paths = readPaths(g, pathsFile, pathIDs);

	switch (opt::mode) {
	  case opt::OVERLAP:
		// Find overlapping paths, do not assemble.
		addPathOverlapEdges(g, paths, pathIDs,
				findOverlaps(g, paths));
		paths.clear();
		if (opt::graphPath.empty())
			opt::graphPath = "-";
		break;

	  case opt::ASSEMBLE:
		// Assemble overlapping paths.
		assembleOverlappingPaths(g, paths, pathIDs);
		break;

	  case opt::TRIM:
		// Trim overlapping paths.
		trimOverlaps(g, paths);
		// Remove paths consisting of a single contig.
		for_each_if(paths.begin(), paths.end(),
				mem_fun_ref(&ContigPath::clear),
				compose1(
					bind2nd(equal_to<ContigPath::size_type>(), 1),
					mem_fun_ref(&ContigPath::size)));
		// Add the paths to the graph.
		addPathOverlapEdges(g, paths, pathIDs, Overlaps());
		break;
	}

	// Output the paths.
	for (Paths::const_iterator it = paths.begin();
			it != paths.end(); ++it) {
		if (it->empty())
			continue;
		assert(it->size() != 1);
		cout << pathIDs[it - paths.begin()] << '\t' << *it << '\n';
	}
	assert(cout.good());

	// Output the graph.
	if (!opt::graphPath.empty()) {
		ofstream fout;
		ostream& out = opt::graphPath == "-" ? cout
		   	: (fout.open(opt::graphPath.c_str()), fout);
		assert_good(out, opt::graphPath);
		write_graph(out, g, PROGRAM, commandLine);
		assert_good(out, opt::graphPath);
	}

	// Output the repeat contigs.
	if (!opt::repeatContigs.empty()) {
		sort(s_trimmedContigs.begin(), s_trimmedContigs.end());
		s_trimmedContigs.erase(
				unique(s_trimmedContigs.begin(),
					s_trimmedContigs.end()), s_trimmedContigs.end());
		ofstream out(opt::repeatContigs.c_str());
		assert_good(out, opt::repeatContigs);
		for (vector<ContigID>::const_iterator it
				= s_trimmedContigs.begin();
				it != s_trimmedContigs.end(); ++it)
			out << ContigID(*it) << '\n';
		assert_good(out, opt::repeatContigs);
	}

	return 0;
}
