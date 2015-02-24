#ifndef CONNECTPAIRS_H
#define CONNECTPAIRS_H

#include "DBGBloomAlgorithms.h"
#include "Bloom/CascadingBloomFilter.h"
#include "DataLayer/FastaInterleave.h"
#include "Graph/BidirectionalBFS.h"
#include "Graph/ConstrainedBidiBFSVisitor.h"
#include "Graph/ExtendPath.h"
#include "Align/alignGlobal.h"
#include "Graph/DefaultColorMap.h"
#include "Graph/DotIO.h"
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <limits>
#include <fstream>

#if _OPENMP
# include <omp.h>
#endif

struct ConnectPairsResult
{
	unsigned k;
	std::string readNamePrefix;
	PathSearchResult pathResult;
	std::vector<FastaRecord> mergedSeqs;
	FastaRecord consensusSeq;
	bool foundStartKmer;
	bool foundGoalKmer;
	unsigned startKmerPos;
	unsigned goalKmerPos;
	unsigned long long numNodesVisited;
	unsigned maxActiveBranches;
	unsigned maxDepthVisitedForward;
	unsigned maxDepthVisitedReverse;
	unsigned pathMismatches;
	unsigned readMismatches;
	size_t memUsage;

	ConnectPairsResult() :
		k(0),
		pathResult(NO_PATH),
		foundStartKmer(false),
		foundGoalKmer(false),
		startKmerPos(NO_MATCH),
		goalKmerPos(NO_MATCH),
		numNodesVisited(0),
		maxActiveBranches(0),
		maxDepthVisitedForward(0),
		maxDepthVisitedReverse(0),
		pathMismatches(0),
		readMismatches(0),
		memUsage(0)
	{}

	static std::ostream& printHeaders(std::ostream& out)
	{
		out << "k\t"
			<< "read_id" << "\t"
			<< "search_result" << "\t"
			<< "num_paths" << "\t"
			<< "path_lengths" << "\t"
			<< "start_kmer_pos" << "\t"
			<< "end_kmer_pos" << "\t"
			<< "nodes_visited" << "\t"
			<< "max_breadth" << "\t"
			<< "max_depth_forward" << "\t"
			<< "max_depth_reverse" << "\t"
			<< "path_mismatches" << "\t"
			<< "read_mismatches" << "\t"
			<< "mem_usage" << "\n";
		return out;
	}

	friend std::ostream& operator <<(std::ostream& out,
		const ConnectPairsResult& o)
	{
		out << o.k << '\t'
			<< o.readNamePrefix << "\t"
			<< PathSearchResultLabel[o.pathResult] << "\t"
			<< o.mergedSeqs.size() << "\t";
		if (o.mergedSeqs.size() == 0) {
			out << "NA" << "\t";
		} else {
			for (unsigned i = 0; i < o.mergedSeqs.size(); i++) {
				out << o.mergedSeqs[i].seq.size();
				if (i < o.mergedSeqs.size() - 1)
					out << ",";
			}
			out << "\t";
		}
		if (o.startKmerPos == NO_MATCH)
			out << "NA\t";
		else
			out << o.startKmerPos << "\t";
		if (o.goalKmerPos == NO_MATCH)
			out << "NA\t";
		else
			out << o.goalKmerPos << "\t";
		out << o.numNodesVisited << "\t"
			<< o.maxActiveBranches << "\t"
			<< o.maxDepthVisitedForward << "\t"
			<< o.maxDepthVisitedReverse << "\t"
			<< o.pathMismatches << "\t"
			<< o.readMismatches << "\t"
			<< o.memUsage << "\n";

		return out;
	}
};

struct ConnectPairsParams {

	unsigned minMergedSeqLen;
	unsigned maxMergedSeqLen;
	unsigned maxPaths;
	unsigned maxBranches;
	unsigned maxPathMismatches;
	unsigned maxReadMismatches;
	unsigned kmerMatchesThreshold;
	bool fixErrors;
	bool maskBases;
	size_t memLimit;
	std::string dotPath;
	std::ofstream* dotStream;

	ConnectPairsParams() :
		minMergedSeqLen(0),
		maxMergedSeqLen(1000),
		maxPaths(NO_LIMIT),
		maxBranches(NO_LIMIT),
		maxPathMismatches(NO_LIMIT),
		maxReadMismatches(NO_LIMIT),
		kmerMatchesThreshold(1),
		fixErrors(false),
		maskBases(false),
		memLimit(std::numeric_limits<std::size_t>::max()),
		dotStream(NULL)
	{}

};

static inline void colorPath(HashGraph<Kmer>& graph, unsigned k,
	const Sequence& seq, const std::string& color,
	bool addEdges = true)
{
	KmerIterator it(seq, k);
	if (it != it.end()) {
		graph.set_vertex_color(*it, color);
		Kmer prev = *it;
		++it;
		for(; it != it.end(); prev=*it, ++it) {
			if (addEdges)
				add_edge(prev, *it, graph);
			graph.set_vertex_color(*it, color);
		}
	}
}

/** Write a color-coded traversal graph to a DOT file. */
static inline void writeDot(
	HashGraph<Kmer>& traversalGraph,
	unsigned k,
	const FastaRecord& read1,
	const FastaRecord& read2,
	const ConnectPairsParams& params,
	const ConnectPairsResult& result)
{
	const std::string pathColor("darkgreen");
	const std::string solutionColor("green");
	const std::string read1Color("blue");
	const std::string read2Color("red");

	// color kmers for the paths / consensus

	const std::vector<FastaRecord>& paths = result.mergedSeqs;

	if (paths.size() == 1) {
		colorPath(traversalGraph, k, paths.front(), solutionColor);
	} else if (paths.size() > 1) {
		for (unsigned i = 0; i < paths.size(); i++)
			colorPath(traversalGraph, k, paths.at(i), pathColor);
		colorPath(traversalGraph, k, result.consensusSeq, solutionColor);
	}

	// color the reads

	colorPath(traversalGraph, k, read1.seq, read1Color);
	colorPath(traversalGraph, k, reverseComplement(read2.seq),
		read2Color);

	// write out the dot file

	// GraphViz utils don't like colons in graph names
	std::string graphName = result.readNamePrefix;
	std::replace(graphName.begin(), graphName.end(), ':', '_');

	write_dot(*params.dotStream, traversalGraph, graphName);
	assert_good(*params.dotStream, params.dotPath);
};

template <typename Graph>
static inline ConnectPairsResult connectPairs(
	unsigned k,
	const FastaRecord& read1,
	const FastaRecord& read2,
	const Graph& g,
	const ConnectPairsParams& params)
{
	ConnectPairsResult result;
	result.k = k;
	result.readNamePrefix = read1.id.substr(0, read1.id.find_last_of("/"));

	if (!isReadNamePair(read1.id, read2.id)) {
#pragma omp critical(cerr)
		std::cerr << "error: name mismatch between paired end reads.\n"
			<< "Read 1: " << read1.id << "\n"
			<< "Read 2: " << read2.id << "\n";
		exit(EXIT_FAILURE);
	}

	if (read1.seq.length() < k || read2.seq.length() < k) {
		result.pathResult = NO_PATH;
		return result;
	}

	unsigned startKmerPos = getStartKmerPos(read1, k, FORWARD, g,
		params.kmerMatchesThreshold);

	unsigned goalKmerPos = getStartKmerPos(read2, k, FORWARD, g,
		params.kmerMatchesThreshold);

	const FastaRecord* pRead1 = &read1;
	const FastaRecord* pRead2 = &read2;
	FastaRecord correctedRead1;
	FastaRecord correctedRead2;
	size_t unused;

	if (startKmerPos == NO_MATCH && params.fixErrors) {
		correctedRead1 = read1;
		if (correctSingleBaseError(g, k, correctedRead1, unused)) {
			startKmerPos = getStartKmerPos(correctedRead1, k, FORWARD, g,
				params.kmerMatchesThreshold);
			assert(startKmerPos != NO_MATCH);
			pRead1 = &correctedRead1;
		}
	}

	if (goalKmerPos == NO_MATCH && params.fixErrors) {
		correctedRead2 = read2;
		if (correctSingleBaseError(g, k, correctedRead2, unused)) {
			goalKmerPos = getStartKmerPos(correctedRead2, k, FORWARD, g,
				params.kmerMatchesThreshold);
			assert(goalKmerPos != NO_MATCH);
			pRead2 = &correctedRead2;
		}
	}

	if (startKmerPos == NO_MATCH || goalKmerPos == NO_MATCH) {
		result.pathResult = NO_PATH;
		return result;
	} else {
		result.startKmerPos = startKmerPos;
		result.foundStartKmer = true;
		result.goalKmerPos = goalKmerPos;
		result.foundGoalKmer = true;
	}

	Kmer startKmer(pRead1->seq.substr(startKmerPos, k));
	Kmer goalKmer(pRead2->seq.substr(goalKmerPos, k));
	goalKmer.reverseComplement();

	unsigned maxPathLen = params.maxMergedSeqLen - k + 1 - startKmerPos - goalKmerPos;
	assert(maxPathLen <= params.maxMergedSeqLen - k + 1);

	unsigned minPathLen = (unsigned)std::max((int)0,
			(int)(params.minMergedSeqLen - k + 1 - startKmerPos - goalKmerPos));
	// do not allow merged seqs that are shorter than the reads
	minPathLen = std::max(minPathLen, (unsigned)std::max(
				pRead1->seq.length() - k + 1 - startKmerPos,
				pRead2->seq.length() - k + 1 - goalKmerPos));

	ConstrainedBidiBFSVisitor<Graph> visitor(g, startKmer, goalKmer,
			params.maxPaths, minPathLen, maxPathLen, params.maxBranches,
			params.memLimit);
	bidirectionalBFS(g, startKmer, goalKmer, visitor);

	std::vector< Path<Kmer> > paths;
	result.pathResult = visitor.pathsToGoal(paths);
	result.numNodesVisited = visitor.getNumNodesVisited();
	result.maxActiveBranches = visitor.getMaxActiveBranches();
	result.maxDepthVisitedForward = visitor.getMaxDepthVisited(FORWARD);
	result.maxDepthVisitedReverse = visitor.getMaxDepthVisited(REVERSE);
	result.memUsage = visitor.approxMemUsage();

	// write traversal graph to dot file (-d option)

	if (result.pathResult == FOUND_PATH) {

		// build sequences for connecting paths

		std::string seqPrefix = pRead1->seq.substr(0, startKmerPos);
		std::string seqSuffix = reverseComplement(pRead2->seq.substr(0, goalKmerPos));
		for (unsigned i = 0; i < paths.size(); i++) {
			FastaRecord mergedSeq;
			std::stringstream index;
			index << i;
			assert(index);
			mergedSeq.id = result.readNamePrefix + "_" + index.str();
			mergedSeq.seq = seqPrefix + pathToSeq(paths[i]) + seqSuffix;
			result.mergedSeqs.push_back(mergedSeq);
		}

		// calc consensus seq and mismatch stats

		if (paths.size() == 1) {

			result.readMismatches =
				maskNew(read1, read2, result.mergedSeqs.front(), params.maskBases);

		} else {

			NWAlignment aln;
			unsigned matches, size;
			boost::tie(matches, size) = align(result.mergedSeqs, aln);
			assert(size >= matches);
			result.pathMismatches = size - matches;

			result.consensusSeq.id = result.readNamePrefix;
			result.consensusSeq.seq = aln.match_align;
			result.readMismatches =
				maskNew(read1, read2, result.consensusSeq, params.maskBases);

		}

	}

	if (!params.dotPath.empty()) {
		HashGraph<Kmer> traversalGraph;
		visitor.getTraversalGraph(traversalGraph);
		writeDot(traversalGraph, k, read1, read2, params, result);
	}

#if 0
# pragma omp critical(cerr)
	std::cerr << result;
#endif

	return result;
}

/**
 * Reason a sequence could not be extended uniquely within
 * the de Bruijn graph; or if the sequence could be extended,
 * the reason we stopped extending.
 */
enum ExtendSeqResult {
	/*
	 * could not find a start kmer in Bloom filter for
	 * path traversal
	 */
	ES_NO_START_KMER,
	/* start kmer had no neighbours */
	ES_DEAD_END,
	/* start kmer had two or more branches */
	ES_BRANCHING_POINT,
	/*
	 * we did not make it from the start kmer to the
	 * beginning/end of the input sequence, because
	 * we hit a dead end.
	 */
	ES_INTERNAL_DEAD_END,
	/*
	 * we did not make it from the start kmer to the
	 * beginning/end of the input sequence, because
	 * we hit a branching point.
	 */
	ES_INTERNAL_BRANCHING_POINT,
	/*
	 * we successfully extended the input sequence
	 * and stopped extending at a branching point.
	 */
	ES_EXTENDED_TO_BRANCHING_POINT,
	/*
	 * we successfully extended the input sequence
	 * and stopped extending at a dead end.
	 */
	ES_EXTENDED_TO_DEAD_END
};

/**
 * Extend a sequence up to the next dead end or branching point in the
 * de Bruijn graph.
 *
 * @param seq sequence to be extended (modified by this function)
 * @param dir direction to extend (FORWARD or REVERSE)
 * @param trimLen ignore branches less than or equal to this length
 * @param k kmer size of de Bruijn graph
 * @param g de Bruijn graph
 * @return ExtendSeqResult (ES_NO_START_KMER, ES_DEAD_END,
 * ES_BRANCHING_POINT, ES_EXTENDED_TO_BRANCHING_POINT,
 * ES_EXTENDED_TO_DEAD_END)
 */
template <typename Graph>
static inline ExtendSeqResult extendSeq(Sequence& seq, Direction dir,
	unsigned k, const Graph& g, unsigned trimLen=0, bool maskNew=true)
{
	if (seq.length() < k)
		return ES_NO_START_KMER;

	size_t origSeqLen = seq.length();

	/* choose a start kmer for the path */

	const unsigned minConsecutiveMatches = 3;
	unsigned startKmerPos = NO_MATCH;
	for (int i = minConsecutiveMatches; i >= 0; i--) {
		startKmerPos = getStartKmerPos2(seq, k, dir, g, i);
		if (startKmerPos != NO_MATCH)
			break;
	}
	if (startKmerPos == NO_MATCH)
		return ES_NO_START_KMER;

	/* initialize the path to be extended */
	std::string kmerStr = seq.substr(startKmerPos, k);
	/* Kmer class doesn't like lowercase chars */
	std::transform(kmerStr.begin(), kmerStr.end(), kmerStr.begin(), ::toupper);
	Kmer startKmer(kmerStr);
	Path<Kmer> path;
	path.push_back(startKmer);

	/*
	 * extend the path up to a dead end or a branching point
	 * within the de Bruijn graph.
	 */

	PathExtensionResult pathResult = extendPath(path, dir, g, trimLen);

	/*
	 * graft path extension onto original input sequence
	 */

	if (pathResult == EXTENDED_TO_DEAD_END ||
		pathResult == EXTENDED_TO_BRANCHING_POINT) {

		/* we expect the start kmer plus some extension */
		assert(path.size() > 1);

		std::string pathSeq = pathToSeq(path);
		std::string::iterator src = pathSeq.begin();
		std::string::iterator dest;

		if (dir == REVERSE) {
			int prefixLen = pathSeq.length() -
				(startKmerPos + k - 1);
			if (prefixLen > 0) {
				/* seq extension goes beyond start of read */
				seq.insert(0, prefixLen, 'N');
			}
			dest = seq.begin();
		} else {
			assert(dir == FORWARD);
			int suffixLen = pathSeq.length() -
				(seq.length() - startKmerPos);
			if (suffixLen > 0) {
				/* seq extension goes beyond end of read */
				seq.insert(seq.length(), suffixLen, 'N');
			}
			dest = seq.begin() + startKmerPos;
		}

		for (; src != pathSeq.end(); ++src, ++dest) {
			assert(dest != seq.end());
			if (maskNew && *src != *dest)
				*src = tolower(*src);
			*dest = *src;
		}

	}

	/* translate and return result code */

	switch (pathResult)
	{
	case EXTENDED_TO_DEAD_END:
		if (seq.length() > origSeqLen)
			return ES_EXTENDED_TO_DEAD_END;
		else
			return ES_INTERNAL_DEAD_END;
	case EXTENDED_TO_BRANCHING_POINT:
		if (seq.length() > origSeqLen)
			return ES_EXTENDED_TO_BRANCHING_POINT;
		else
			return ES_INTERNAL_BRANCHING_POINT;
	case DEAD_END:
		return ES_DEAD_END;
	case BRANCHING_POINT:
		return ES_BRANCHING_POINT;
	default:
		break;
	}

	/* should never reach here */
	assert(false);
}

#endif
