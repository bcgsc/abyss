/**
 * Algorithms for a de Bruijn Graph using a Bloom filter
 * Copyright 2013 Canada's Michael Smith Genome Science Centre
 */
#ifndef DBGBLOOMALGORITHMS_H
#define DBGBLOOMALGORITHMS_H 1

#include "Common/Kmer.h"
#include "Common/Warnings.h"
#include "DBGBloom.h"
#include "Common/StringUtil.h"
#include "Common/Sequence.h"
#include "DataLayer/FastaReader.h"
#include "Graph/DefaultColorMap.h"
#include "Graph/Path.h"
#include "Graph/BidirectionalBFS.h"
#include "Graph/ConstrainedBidiBFSVisitor.h"
#include <climits>
#include <algorithm> // for std::max

#if _OPENMP
# include <omp.h>
#endif

#define NO_MATCH UINT_MAX

struct SearchResult
{
	std::string readNamePrefix;
	PathSearchResult pathResult;
	std::vector<FastaRecord> mergedSeqs;
	bool foundStartKmer;
	bool foundGoalKmer;
	unsigned startKmerPos;
	unsigned goalKmerPos;
	unsigned long long numNodesVisited;
	unsigned maxActiveBranches;
	unsigned maxDepthVisitedForward;
	unsigned maxDepthVisitedReverse;

	SearchResult() :
		pathResult(NO_PATH),
		foundStartKmer(false),
		foundGoalKmer(false),
		startKmerPos(NO_MATCH),
		goalKmerPos(NO_MATCH),
		numNodesVisited(0),
		maxActiveBranches(0),
		maxDepthVisitedForward(0),
		maxDepthVisitedReverse(0) {}

	static std::ostream& printHeaders(std::ostream& out)
	{
		out << "read_id" << "\t"
			<< "search_result" << "\t"
			<< "num_paths" << "\t"
			<< "path_lengths" << "\t"
			<< "start_kmer_pos" << "\t"
			<< "end_kmer_pos" << "\t"
			<< "nodes_visited" << "\t"
			<< "max_breadth" << "\t"
			<< "max_depth_forward" << "\t"
			<< "max_depth_reverse" << "\n";
		return out;
	}

	friend std::ostream& operator <<(std::ostream& out,
		const SearchResult& o)
	{
		out << o.readNamePrefix << "\t"
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
		out << o.startKmerPos << "\t"
			<< o.goalKmerPos << "\t"
			<< o.numNodesVisited << "\t"
			<< o.maxActiveBranches << "\t"
			<< o.maxDepthVisitedForward << "\t"
			<< o.maxDepthVisitedReverse << "\n";
		return out;
	}
};

static inline Sequence pathToSeq(Path<Kmer> path)
{
	Sequence seq;
	assert(path.size() > 0);
	seq.append(path[0].str());
	for (unsigned i = 1; i < path.size(); i++)
		seq.append(1, path[i].getLastBaseChar());
	return seq;
}

static inline unsigned getStartKmerPos(unsigned k,
	const FastaRecord& read, const DBGBloom& g,
	bool rc = false)
{
	// build a vector indicating whether each kmer is a match
	// The vector intentionally has an extra false element at
	// the end, for the second loop below.

	const std::string& seq = read.seq;
	std::vector<bool> match(seq.length() - k + 2, false);
	bool foundMatch = false;
	for (unsigned i = 0; i < seq.length() - k + 1; i++) {
		std::string kmerStr = seq.substr(i,k);
		if (kmerStr.find_first_not_of("AGCTagct") != std::string::npos)
			continue;
		Kmer kmer(kmerStr);
		if (rc)
			kmer.reverseComplement();
		if (graph_traits<DBGBloom>::vertex_exists(kmer, g)) {
			foundMatch = true;
			match[i] = true;
		}
	}
	if (!foundMatch)
		return NO_MATCH;

	// find the longest string of matches

	unsigned maxMatchLength = 0;
	unsigned maxMatchPos = 0;
	unsigned matchLength = 0;
	unsigned matchPos = 0;
	bool matchPosSet = false;
	for (unsigned i = 0; i < match.size(); i++) {
		if (match[i]) {
			if (!matchPosSet) {
				matchPos = i;
				matchPosSet = true;
			}
			matchLength++;
		} else {
			// Note: match has an extra false element at the end,
			// so this else block will get executed at least once.
			if (matchLength >= maxMatchLength) {
				maxMatchPos = matchPos;
				maxMatchLength = matchLength;
			}
			matchLength = 0;
			matchPosSet = false;
		}
	}
	assert(maxMatchLength > 0);

	// return the kmer closest to the gap between read pairs
	return maxMatchPos + maxMatchLength - 1;
}

static inline SearchResult connectPairs(
	unsigned k,
	const FastaRecord& read1,
	const FastaRecord& read2,
	const DBGBloom& g,
	unsigned maxPaths = 2,
	unsigned minMergedSeqLen = 0,
	unsigned maxMergedSeqLen = NO_LIMIT,
	unsigned maxBranches = NO_LIMIT)
{
	SearchResult result;

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

	unsigned startKmerPos = getStartKmerPos(k, read1, g, false);
	unsigned goalKmerPos = getStartKmerPos(k, read2, g, true);

	if (startKmerPos != NO_MATCH) {
		result.startKmerPos = startKmerPos;
		result.foundStartKmer = true;
	}

	if (goalKmerPos != NO_MATCH) {
		result.goalKmerPos = goalKmerPos;
		result.foundGoalKmer = true;
	}

	if (startKmerPos == NO_MATCH || goalKmerPos == NO_MATCH) {
		result.pathResult = NO_PATH;
		return result;
	}

	Kmer startKmer(read1.seq.substr(startKmerPos, k));
	Kmer goalKmer(read2.seq.substr(goalKmerPos, k));
	goalKmer.reverseComplement();

	unsigned maxPathLen = maxMergedSeqLen - k + 1 - startKmerPos - goalKmerPos;
	assert(maxPathLen <= maxMergedSeqLen - k + 1);
	unsigned longestReadLen = std::max(read1.seq.length(), read2.seq.length());
	// note: do not allow paths that are shorter than the reads
	unsigned minPathLen = (unsigned)std::max((int)(longestReadLen - k + 1),
			(int)(minMergedSeqLen - k + 1 - startKmerPos - goalKmerPos));

	ConstrainedBidiBFSVisitor<DBGBloom> visitor(g, startKmer, goalKmer, maxPaths,
			minPathLen, maxPathLen, maxBranches);
	bidirectionalBFS(g, startKmer, goalKmer, visitor);

	std::vector< Path<Kmer> > paths;
	result.readNamePrefix = read1.id.substr(0, read1.id.find_last_of("/"));
	result.pathResult = visitor.pathsToGoal(paths);
	result.numNodesVisited = visitor.getNumNodesVisited();
	result.maxActiveBranches = visitor.getMaxActiveBranches();
	result.maxDepthVisitedForward = visitor.getMaxDepthVisited(FORWARD);
	result.maxDepthVisitedReverse = visitor.getMaxDepthVisited(REVERSE);

	if (result.pathResult == FOUND_PATH) {
		std::string seqPrefix = read1.seq.substr(0, startKmerPos);
		std::string seqSuffix = reverseComplement(read2.seq.substr(0, goalKmerPos));
		for (unsigned i = 0; i < paths.size(); i++) {
			FastaRecord mergedSeq;
			mergedSeq.id = result.readNamePrefix;
			mergedSeq.seq = seqPrefix + pathToSeq(paths[i]) + seqSuffix;
			result.mergedSeqs.push_back(mergedSeq);
		}
	}

#if 0
# pragma omp critical(cerr)
	std::cerr << result;
#endif

	return result;
}

#endif
