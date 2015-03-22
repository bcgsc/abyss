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
#include "Common/Sequence.h"
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

static inline unsigned getHeadKmerPos(const Sequence& seq, Direction dir,
	unsigned k)
{
	return (dir == FORWARD) ?  seq.length() - k : 0;
}

static inline Kmer getHeadKmer(const Sequence& seq, Direction dir,
	unsigned k)
{
	return Kmer(seq.substr(getHeadKmerPos(seq, dir, k), k));
}

template <typename Graph>
static inline bool extendSeqThroughBubble(Sequence& seq,
	Direction dir, unsigned startKmerPos, unsigned k,
	const Graph& g, unsigned trimLen=0, bool maskNew=false)
{
	assert(seq.length() >= k);
	assert(dir == FORWARD || dir == REVERSE);

	/*
	 * unhandled case: bubble is contained entirely
	 * within input sequence.
	 */
	unsigned bubbleSeqLen = 2*k + 1;

	if (dir == FORWARD &&
		startKmerPos + bubbleSeqLen <= seq.length()) {
		return false;
	} else if (dir == REVERSE &&
		bubbleSeqLen <= startKmerPos) {
		return false;
	}

	Kmer head(seq.substr(startKmerPos, k));
	std::vector<Kmer> buds = trueBranches(head, dir, g, trimLen);

	/* more than two branches -- not a simple bubble */
	if (buds.size() != 2)
		return false;

	Path<Kmer> path1, path2;
	path1.push_back(buds.front());
	path2.push_back(buds.back());
	extendPath(path1, dir, g, trimLen, k+2);
	extendPath(path2, dir, g, trimLen, k+2);

	/* paths lengths not k+1 -- not a simple bubble */
	if (path1.size() != k+2 || path2.size() != k+2)
		return false;

	Kmer head1, head2;
	if (dir == FORWARD) {
		head1 = path1.back();
		head2 = path2.back();
	} else {
		assert(dir == REVERSE);
		head1 = path1.front();
		head2 = path2.front();
	}

	/* paths don't reconnect -- not a simple bubble */
	if (head1 != head2)
		return false;

	NWAlignment alignment;
	alignPair(pathToSeq(path1), pathToSeq(path2), alignment);
	Sequence& consensus = alignment.match_align;

	if (dir == FORWARD) {
		overlaySeq(consensus, seq, seq.length()-k, maskNew);
	} else {
		overlaySeq(consensus, seq, -consensus.length()+k, maskNew);
	}

	return true;
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
	ES_NO_START_KMER=0,
	/* start kmer had no neighbours */
	ES_DEAD_END,
	/* start kmer had two or more branches */
	ES_BRANCHING_POINT,
	/* start kmer was part of a cycle */
	ES_CYCLE,
	/* input seq was already max length or less */
	ES_LENGTH_LIMIT,
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
	 * we did not make from the start kmer to the
	 * beginning/end of the input sequence, because
	 * we hit a cycle.
	 */
	ES_INTERNAL_CYCLE,
	/*
	 * we successfully extended the input sequence
	 * and stopped extending at a dead end.
	 */
	ES_EXTENDED_TO_DEAD_END,
	/*
	 * we successfully extended the input sequence
	 * and stopped extending at a branching point.
	 */
	ES_EXTENDED_TO_BRANCHING_POINT,
	/*
	 * we successfully extended the input sequence
	 * and stopped when we hit a cycle.
	 */
	ES_EXTENDED_TO_CYCLE,
	/*
	 * we successfully extended the input sequence
	 * the given length limit
	 */
	ES_EXTENDED_TO_LENGTH_LIMIT
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
	unsigned startKmerPos, unsigned k, const Graph& g,
	unsigned maxLen=NO_LIMIT, unsigned trimLen=0,
	bool maskNew=false)
{
	if (seq.length() < k)
		return ES_NO_START_KMER;

	assert(startKmerPos < seq.length()-k+1);

	if (maxLen < seq.length())
		maxLen = seq.length();

	size_t origSeqLen = seq.length();

#if 0
Sequence origSeq = seq;
#endif

	/* initialize the path to be extended */
	std::string kmerStr = seq.substr(startKmerPos, k);
	if (kmerStr.find_first_not_of("AGCTagct") !=
		std::string::npos)
		return ES_NO_START_KMER;

	Kmer startKmer(kmerStr);
	Path<Kmer> path;
	path.push_back(startKmer);
	PathExtensionResult pathResult = LENGTH_LIMIT;

	bool done = false;
	while (!done && seq.length() < maxLen) {

		/*
		 * extend the path up to a dead end or branching point
		 * in the de Bruijn graph.
		 */
		unsigned maxPathLen;
		if (maxLen == NO_LIMIT) {
			maxPathLen = NO_LIMIT;
		} else if (dir == FORWARD) {
			maxPathLen = std::max((int)1,
				(int)(maxLen - startKmerPos - k + 1));
		} else {
			assert(dir == REVERSE);
			maxPathLen = std::max((int)1,
				(int)(maxLen - seq.length() + startKmerPos + 1));
		}

		pathResult = extendPath(path, dir, g, trimLen,
			maxPathLen);

		/*
		 * graft path extension onto original input sequence
		 */

		size_t lengthBeforeOverlay = seq.length();
		if (pathResult == EXTENDED_TO_DEAD_END ||
			pathResult == EXTENDED_TO_BRANCHING_POINT ||
			pathResult == EXTENDED_TO_CYCLE ||
			pathResult == EXTENDED_TO_LENGTH_LIMIT)
		{
			std::string pathSeq = pathToSeq(path);
			if (dir == FORWARD) {
				overlaySeq(pathSeq, seq, startKmerPos, maskNew);
			} else {
				assert(dir == REVERSE);
				overlaySeq(pathSeq, seq, -pathSeq.length() + startKmerPos + k,
					maskNew);
			}
		}
		size_t lengthAfterOverlay = seq.length();

		/*
		 * extend through simple bubbles
		 */
		done = true;
		if (lengthAfterOverlay < maxLen &&
			(pathResult == BRANCHING_POINT ||
			pathResult == EXTENDED_TO_BRANCHING_POINT)) {
			if (dir == FORWARD) {
				startKmerPos = startKmerPos + path.size() - 1;
				assert(startKmerPos < seq.length() - k + 1);
			} else {
				assert(dir == REVERSE);
				if (lengthAfterOverlay > lengthBeforeOverlay) {
					startKmerPos = 0;
				} else {
					assert(startKmerPos >= path.size() - 1);
					startKmerPos = startKmerPos - path.size() + 1;
				}
			}
			if (extendSeqThroughBubble(seq, dir, startKmerPos,
				k, g, trimLen, maskNew)) {
				if (seq.length() > maxLen) {
					if (dir == FORWARD) {
						seq = seq.substr(0, maxLen);
					} else {
						assert(dir == REVERSE);
						seq = seq.substr(seq.length() - maxLen,
							maxLen);
					}
				} else {
					path.clear();
					path.push_back(getHeadKmer(seq, dir, k));
					startKmerPos = getHeadKmerPos(seq, dir, k);
					done = false;
				}
			}
		}

	} /* while (!done && seq.length() < maxLen) */

	/* translate and return result code */

	ExtendSeqResult result;

	switch (pathResult)
	{
	case EXTENDED_TO_DEAD_END:
		if (seq.length() > origSeqLen)
			result = ES_EXTENDED_TO_DEAD_END;
		else
			result = ES_INTERNAL_DEAD_END;
		break;
	case EXTENDED_TO_BRANCHING_POINT:
		if (seq.length() > origSeqLen)
			result = ES_EXTENDED_TO_BRANCHING_POINT;
		else
			result = ES_INTERNAL_BRANCHING_POINT;
		break;
	case EXTENDED_TO_CYCLE:
		if (seq.length() > origSeqLen)
			result = ES_EXTENDED_TO_CYCLE;
		else
			result = ES_INTERNAL_CYCLE;
	case DEAD_END:
		result = ES_DEAD_END;
		break;
	case BRANCHING_POINT:
		result = ES_BRANCHING_POINT;
		break;
	case CYCLE:
		result = ES_CYCLE;
		break;
	case LENGTH_LIMIT:
		result = ES_LENGTH_LIMIT;
		break;
	case EXTENDED_TO_LENGTH_LIMIT:
		result = ES_EXTENDED_TO_LENGTH_LIMIT;
		break;
	default:
		/* all other cases should be handled above */
		assert(false);
	}

#if 0
#pragma omp critical(cerr)
{
	std::cerr << "-----\n";
	std::cerr << "extend dir: " << ((dir == FORWARD) ? "forward" : "reverse") << "\n";
	std::cerr << "startKmerPos: " << startKmerPos << "\n";
	std::cerr << "result: " << result << "\n";
	std::cerr << "orig len: " << origSeq.length() << "\n";
	std::cerr << "extended len: " << seq.length() << "\n";
	if (result == ES_EXTENDED_TO_DEAD_END ||
		result == ES_EXTENDED_TO_BRANCHING_POINT) {
		assert(seq.length() > origSeq.length());
		std::string pathSeq = pathToSeq(path);
		std::cerr << "origSeq:     ";
		unsigned padding = seq.length() - origSeq.length();
		if (dir == REVERSE) {
			for (unsigned i = 0; i < padding; ++i)
				std::cerr << "-";
		}
		std::cerr << origSeq;
		if (dir == FORWARD) {
			for (unsigned i = 0; i < padding; ++i)
				std::cerr << "-";
		}
		std::cerr << "\n";
		std::cerr << "pathSeq:     ";
		if (dir == FORWARD) {
			for (unsigned i = 0; i < startKmerPos; ++i)
				std::cerr << " ";
		}
		std::cerr << pathSeq << "\n";
		std::cerr << "extendedSeq: " << seq << "\n";
	}
	std::cerr << "-----\n";
}
#endif

	return result;
}

/**
 * Correct the given sequence using the Bloom filter de Bruijn
 * graph.  The correction is performed by finding the longest
 * stretch of good kmers in the sequence and extending that
 * region both left and right.
 */
template <typename Graph>
static inline bool correctAndExtendSeq(Sequence& seq,
	unsigned k, const Graph& g, unsigned maxLen=NO_LIMIT,
	unsigned trimLen=0, bool maskNew=false)
{
	if (seq.size() < k)
		return false;

	if (maxLen < seq.length())
		maxLen = seq.length();

	/*
	 * find longest stretch of contiguous kmers
	 * in de Bruijn graph
	 */

	const unsigned UNSET = UINT_MAX;
	unsigned matchStart = UNSET;
	unsigned matchLen = 0;
	unsigned maxMatchLen = 0;
	unsigned maxMatchStart = UNSET;

	for (unsigned i = 0; i < seq.length() - k + 1; ++i) {
		std::string kmerStr = seq.substr(i, k);
		size_t pos = kmerStr.find_first_not_of("AGCTagct");
		if (pos != std::string::npos ||
			!vertex_exists(Kmer(kmerStr), g)) {
			if (matchStart != UNSET &&
				matchLen > maxMatchLen) {
				maxMatchLen = matchLen;
				maxMatchStart = matchStart;
			}
			matchStart = UNSET;
			matchLen = 0;
			if (pos != std::string::npos)
				i += pos;
		} else {
			if (matchStart == UNSET)
				matchStart = i;
			matchLen++;
		}
	}
	if (matchStart != UNSET && matchLen > maxMatchLen) {
		maxMatchStart = matchStart;
		maxMatchLen = matchLen;
	}
	if (maxMatchLen == 0)
		return false;
	assert(maxMatchStart != UNSET);
	assert(maxMatchLen > 0);

	unsigned maxMatchSeqLen = maxMatchLen+k-1;
	unsigned seedSeqLen = std::min(2*k-1, maxMatchSeqLen);

	Sequence correctedSeq = seq.substr(
		maxMatchStart + maxMatchSeqLen - seedSeqLen,
		std::string::npos);

	extendSeq(correctedSeq, REVERSE, correctedSeq.length()-k, k, g, 2*k,
		trimLen, maskNew);
	if (correctedSeq.length() < 2*k)
		return false;

	correctedSeq = correctedSeq.substr(0, k);
	extendSeq(correctedSeq, FORWARD, 0, k, g, 2*k+1, trimLen, maskNew);

	seq = correctedSeq;
	return true;
}

#endif
