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
#include "Common/KmerSet.h"
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
	unsigned searchCost;
	/** alternate connecting sequence(s) for read pair */
	std::vector<Sequence> connectingSeqs;
	/** read pairs joined with alternate connecting sequence(s) */
	std::vector<FastaRecord> mergedSeqs;
	/** consensus sequence for alternate connecting sequences */
	Sequence consensusConnectingSeq;
	/**
	 * consensus sequence for read pairs joined by
	 * alternate connecting sequences
	 */
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
	float pathIdentity;
	unsigned readMismatches;
	float readIdentity;
	size_t memUsage;

	ConnectPairsResult() :
		k(0),
		pathResult(NO_PATH),
		searchCost(0),
		foundStartKmer(false),
		foundGoalKmer(false),
		startKmerPos(NO_MATCH),
		goalKmerPos(NO_MATCH),
		numNodesVisited(0),
		maxActiveBranches(0),
		maxDepthVisitedForward(0),
		maxDepthVisitedReverse(0),
		pathMismatches(0),
		pathIdentity(0.0f),
		readMismatches(0),
		readIdentity(0.0f),
		memUsage(0)
	{}

	static std::ostream& printHeaders(std::ostream& out)
	{
		out << "k\t"
			<< "read_id" << "\t"
			<< "search_result" << "\t"
			<< "search_cost" << "\t"
			<< "num_paths" << "\t"
			<< "path_lengths" << "\t"
			<< "start_kmer_pos" << "\t"
			<< "end_kmer_pos" << "\t"
			<< "nodes_visited" << "\t"
			<< "max_breadth" << "\t"
			<< "max_depth_forward" << "\t"
			<< "max_depth_reverse" << "\t"
			<< "path_mismatches" << "\t"
			<< "path_identity" << "\t"
			<< "read_mismatches" << "\t"
			<< "read_identity" << "\t"
			<< "mem_usage" << "\n";
		return out;
	}

	friend std::ostream& operator <<(std::ostream& out,
		const ConnectPairsResult& o)
	{
		out << o.k << '\t'
			<< o.readNamePrefix << "\t"
			<< PathSearchResultLabel[o.pathResult] << "\t"
			<< o.searchCost << "\t"
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
			<< std::setprecision(3) << o.pathIdentity << "\t"
			<< o.readMismatches << "\t"
			<< std::setprecision(3) << o.readIdentity << "\t"
			<< o.memUsage << "\n";

		return out;
	}
};

struct ConnectPairsParams {

	unsigned minMergedSeqLen;
	unsigned maxMergedSeqLen;
	unsigned maxPaths;
	unsigned maxBranches;
	unsigned maxCost;
	unsigned maxPathMismatches;
	float minPathIdentity;
	unsigned maxReadMismatches;
	float minReadIdentity;
	unsigned kmerMatchesThreshold;
	bool fixErrors;
	bool maskBases;
	bool preserveReads;
	size_t memLimit;
	std::string dotPath;
	std::ofstream* dotStream;

	ConnectPairsParams() :
		minMergedSeqLen(0),
		maxMergedSeqLen(1000),
		maxPaths(NO_LIMIT),
		maxBranches(NO_LIMIT),
		maxCost(NO_LIMIT),
		maxPathMismatches(NO_LIMIT),
		minPathIdentity(0.0f),
		maxReadMismatches(NO_LIMIT),
		minReadIdentity(0.0f),
		kmerMatchesThreshold(1),
		fixErrors(false),
		maskBases(false),
		preserveReads(false),
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

	const unsigned numMatchesThreshold = 3;

	unsigned startKmerPos = getStartKmerPos(read1, k, FORWARD, g,
		numMatchesThreshold, params.preserveReads);

	unsigned goalKmerPos = getStartKmerPos(read2, k, FORWARD, g,
		numMatchesThreshold, params.preserveReads);

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
			params.maxCost, params.memLimit);
	bidirectionalBFS(g, startKmer, goalKmer, visitor);

	std::vector< Path<Kmer> > paths;
	result.pathResult = visitor.pathsToGoal(paths);
	result.searchCost = visitor.getSearchCost();
	result.numNodesVisited = visitor.getNumNodesVisited();
	result.maxActiveBranches = visitor.getMaxActiveBranches();
	result.maxDepthVisitedForward = visitor.getMaxDepthVisited(FORWARD);
	result.maxDepthVisitedReverse = visitor.getMaxDepthVisited(REVERSE);
	result.memUsage = visitor.approxMemUsage();

	if (result.pathResult == FOUND_PATH) {

		/* build sequences for connecting paths */

		std::string seqPrefix, seqSuffix;

		if (params.preserveReads) {
			seqPrefix = pRead1->seq;
			seqSuffix = reverseComplement(pRead2->seq);
			unsigned trimLeft = pRead1->seq.length() - startKmerPos;
			unsigned trimRight = pRead2->seq.length() - goalKmerPos;
			for (unsigned i = 0; i < paths.size(); i++) {
				Sequence connectingSeq = pathToSeq(paths[i]);
				/*
				 * If the input reads overlap, we must fail because
				 * there's no way to preserve the original read
				 * sequences in the merged read (the reads may disagree
				 * in the region of overlap)
				 */
				if (trimLeft + trimRight > connectingSeq.length()) {
					result.pathResult = NO_PATH;
					return result;
				}
				connectingSeq = connectingSeq.substr(trimLeft,
					connectingSeq.length() - trimLeft - trimRight);
				result.connectingSeqs.push_back(connectingSeq);
			}
		} else {
			seqPrefix = pRead1->seq.substr(0, startKmerPos);
			seqSuffix = reverseComplement(pRead2->seq.substr(0, goalKmerPos));
			for (unsigned i = 0; i < paths.size(); i++)
				result.connectingSeqs.push_back(pathToSeq(paths[i]));
		}

		unsigned readPairLength = read1.seq.length() + read2.seq.length();

		if (paths.size() == 1) {

			/* found a unique path between the reads */

			FastaRecord mergedSeq;
			mergedSeq.id = result.readNamePrefix;
			mergedSeq.seq = seqPrefix + result.connectingSeqs.front() + seqSuffix;
			result.readMismatches =
				maskNew(read1, read2, mergedSeq, params.maskBases);
			result.pathIdentity = 100.0f;
			result.readIdentity = 100.0f * (float)(readPairLength -
				result.readMismatches) / readPairLength;

			result.mergedSeqs.push_back(mergedSeq);
			result.consensusSeq = mergedSeq;
			result.consensusConnectingSeq = result.connectingSeqs.front();

		} else {

			/*
			 * multiple paths were found, so build a consensus
			 * sequence using multiple sequence alignment.
			 */

			NWAlignment aln;
			unsigned matches, size;
			boost::tie(matches, size) = align(result.connectingSeqs, aln);
			assert(size >= matches);
			result.pathMismatches = size - matches;
			result.consensusConnectingSeq = aln.match_align;
			result.pathIdentity = 100.0f *
				(float)(result.consensusConnectingSeq.length()
				- result.pathMismatches) / result.consensusConnectingSeq.length();
			result.consensusSeq.id = result.readNamePrefix;
			result.consensusSeq.seq = seqPrefix + result.consensusConnectingSeq +
				seqSuffix;
			result.readMismatches =
				maskNew(read1, read2, result.consensusSeq, params.maskBases);
			result.readIdentity = 100.0f * (float)(readPairLength -
				result.readMismatches) / readPairLength;

			unsigned i = 1;
			for (std::vector<Sequence>::iterator it = result.connectingSeqs.begin();
				it != result.connectingSeqs.end(); ++it) {
				FastaRecord mergedSeq;
				std::ostringstream id;
				id << result.readNamePrefix << '_' << i++;
				mergedSeq.id = id.str();
				mergedSeq.seq = seqPrefix + *it + seqSuffix;
				result.mergedSeqs.push_back(mergedSeq);
			}

		}

		assert(result.connectingSeqs.size() == result.mergedSeqs.size());
	}

	/* write traversal graph to dot file (-d option) */

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
	const Graph& g, unsigned trimLen=0, bool maskNew=false,
	bool preserveSeq=false)
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

	std::string headKmer = seq.substr(startKmerPos, k);
	if (headKmer.find_first_not_of("AGCTagct") !=
		std::string::npos)
		return false;

	/* longest branch of Bloom filter false positives */
	const unsigned fpTrim = 5;

	Kmer head(seq.substr(startKmerPos, k));
	std::vector<Kmer> buds = trueBranches(head, dir, g, trimLen, fpTrim);

	/* more than two branches -- not a simple bubble */
	if (buds.size() != 2)
		return false;

	Path<Kmer> path1, path2;
	if (dir == FORWARD) {
		path1.push_back(head);
		path2.push_back(head);
	}
	path1.push_back(buds.front());
	path2.push_back(buds.back());
	if (dir == REVERSE) {
		path1.push_back(head);
		path2.push_back(head);
	}

	ExtendPathParams params;
	params.trimLen = trimLen;
	params.fpTrim = 5;
	params.maxLen = k + 2;
	params.lookBehind = true;

	extendPath(path1, dir, g, params);
	extendPath(path2, dir, g, params);

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
		if (preserveSeq) {
			/*
			 * make sure bubble extends beyond end of
			 * original sequence
			 */
			assert(startKmerPos + consensus.length()
				> seq.length());
			overlaySeq(consensus.substr(seq.length() - startKmerPos),
				seq, seq.length(), maskNew);
		} else {
			overlaySeq(consensus, seq, startKmerPos, maskNew);
		}
	} else {
		if (preserveSeq) {
			/*
			 * make sure bubble extends beyond end of
			 * original sequence
			 */
			assert(consensus.length() > startKmerPos + k);
			consensus = consensus.substr(0,
				consensus.length() - startKmerPos - k);
			overlaySeq(consensus, seq, -consensus.length(),
				maskNew);
		} else {
			overlaySeq(consensus, seq,
				-consensus.length() + startKmerPos + k, maskNew);
		}
	}

	return true;
}

Path<Kmer> seqToPath(const Sequence& seq, unsigned k)
{
	assert(seq.length() >= k);
	Path<Kmer> path;
	Sequence seqCopy = seq;
	flattenAmbiguityCodes(seqCopy);
	for (unsigned i = 0; i < seq.length() - k + 1; ++i) {
		std::string kmerStr = seq.substr(i, k);
		path.push_back(Kmer(kmerStr));
	}
	return path;
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
	bool maskNew=false, bool popBubbles=true,
	bool preserveSeq=false)
{
	if (seq.length() < k)
		return ES_NO_START_KMER;

	assert(startKmerPos < seq.length() - k + 1);

	if (maxLen < seq.length())
		maxLen = seq.length();

	size_t origSeqLen = seq.length();

	/*
	 * temporarily switch orientation so that REVERSE and FORWARD cases
	 * can be handled in the same way.
	 */

	if (dir == REVERSE) {
		startKmerPos = seq.length() - startKmerPos - k;
		assert(startKmerPos < seq.length() - k + 1);
		seq = reverseComplement(seq);
	}

	/* initialize the path to be extended */

	std::string kmerStr = seq.substr(startKmerPos, k);
	if (kmerStr.find_first_not_of("AGCTagct") !=
		std::string::npos)
		return ES_NO_START_KMER;

	Kmer startKmer(kmerStr);
	Path<Kmer> path;
	path.push_back(startKmer);
	PathExtensionResult pathResult = std::make_pair(0, ER_DEAD_END);

	/* track visited kmers to avoid traversing cycles in an infinite loop */

	KmerSet visited(k);
	visited.loadSeq(seq);

	/* extend through unambiguous paths and simple bubbles */

	bool done = false;
	while (!done && seq.length() < maxLen) {

		/*
		 * extend the path up to the next dead end or branching point
		 * in the de Bruijn graph.
		 */
		unsigned maxPathLen;
		if (maxLen == NO_LIMIT) {
			maxPathLen = NO_LIMIT;
		} else {
			maxPathLen = (unsigned)std::max((int)1,
				(int)(maxLen - startKmerPos - k + 1));
		}

		ExtendPathParams params;
		params.trimLen = trimLen;
		params.fpTrim = 5;
		params.maxLen = maxPathLen;
		params.lookBehind = false;

		pathResult = extendPath(path, FORWARD, g, params);

		/*
		 * give up if we don't at extend beyond end
		 * of existing sequence
		 */
		unsigned overlappingKmers = seq.length() - startKmerPos - k + 1;
		if (path.size() <= overlappingKmers) {
			done = true;
			break;
		}

		/* check for cycle */
		path.erase(path.begin(), path.begin() + overlappingKmers);
		for (Path<Kmer>::iterator it = path.begin();
			it != path.end(); ++it) {
			if (visited.containsKmer(*it)) {
				pathResult.second = ER_CYCLE;
				path.erase(it, path.end());
				break;
			}
			visited.addKmer(*it);
		}

		/*
		 * graft path extension onto original input sequence
		 */
		if (path.size() > 0 && pathResult.first > 0)
		{
			std::string pathSeq = pathToSeq(path);
			if (preserveSeq)
				overlaySeq(pathSeq.substr(k), seq,
					seq.length(), maskNew);
			else
				overlaySeq(pathSeq, seq,
					seq.length() - k + 1, maskNew);
		}

		/*
		 * extend through simple bubbles
		 */
		done = true;
		if (popBubbles && seq.length() < maxLen
			&& pathResult.second == ER_AMBI_OUT) {
			startKmerPos = startKmerPos + path.size() - 1;
			assert(startKmerPos < seq.length() - k + 1);
			if (extendSeqThroughBubble(seq, FORWARD, startKmerPos,
				k, g, trimLen, maskNew, preserveSeq)) {

				/* make sure we don't exceed extension limit */
				if (seq.length() > maxLen)
					seq = seq.substr(0, maxLen);

				/* check for cycle */
				for (unsigned i = startKmerPos + 1;
					i < seq.length() - k + 1; ++i) {
					std::string kmerStr = seq.substr(i, k);
					size_t pos = kmerStr.find_first_not_of("AGCTagct");
					if (pos != std::string::npos) {
						i += pos;
						continue;
					}
					Kmer kmer(kmerStr);
					if (visited.containsKmer(kmer)) {
						pathResult.second = ER_CYCLE;
						seq.erase(i);
						break;
					}
					visited.addKmer(kmer);
				}

				/* set up for another round of extension */
				if (pathResult.second != ER_CYCLE && seq.length() < maxLen) {
					done = false;
					startKmerPos = seq.length() - k;
					path.clear();
					path.push_back(Kmer(seq.substr(startKmerPos)));
				}
			}
		}

	} /* while (!done && seq.length() < maxLen) */

	/* translate and return result code */

	ExtendSeqResult result;

	switch (pathResult.second)
	{
	case ER_DEAD_END:
		if (seq.length() == k)
			result = ES_DEAD_END;
		else if (seq.length() > origSeqLen)
			result = ES_EXTENDED_TO_DEAD_END;
		else
			result = ES_INTERNAL_DEAD_END;
		break;
	case ER_AMBI_IN:
	case ER_AMBI_OUT:
		if (seq.length() == k)
			result = ES_BRANCHING_POINT;
		else if (seq.length() > origSeqLen)
			result = ES_EXTENDED_TO_BRANCHING_POINT;
		else
			result = ES_INTERNAL_BRANCHING_POINT;
		break;
	case ER_CYCLE:
		if (seq.length() == k)
			result = ES_CYCLE;
		else if (seq.length() > origSeqLen)
			result = ES_EXTENDED_TO_CYCLE;
		else
			result = ES_INTERNAL_CYCLE;
		break;
	case ER_LENGTH_LIMIT:
		if (seq.length() == origSeqLen) {
			result = ES_LENGTH_LIMIT;
		} else  {
			assert(seq.length() > origSeqLen);
			result = ES_EXTENDED_TO_LENGTH_LIMIT;
		}
		break;
	default:
		/* all other cases should be handled above */
		assert(false);
	}

	/* switch back to original orientation */
	if (dir == REVERSE)
		seq = reverseComplement(seq);

	return result;
}

template <typename Graph>
static inline bool trimRead(FastqRecord& read,
	unsigned k, const Graph& g)
{
	Sequence& seq = read.seq;

	if (seq.size() < k)
		return false;

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

	read.seq = read.seq.substr(maxMatchStart, maxMatchLen + k - 1);
	if (!read.qual.empty())
		read.qual = read.qual.substr(maxMatchStart, maxMatchLen + k - 1);
	return true;
}

#endif
