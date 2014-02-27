#ifndef CONNECTPAIRS_H
#define CONNECTPAIRS_H

#include "connectpairs/DBGBloomAlgorithms.h"
#include "DataLayer/FastaInterleave.h"
#include "Graph/BidirectionalBFS.h"
#include "Graph/ConstrainedBidiBFSVisitor.h"
#include "Align/alignGlobal.h"
#include "Graph/DefaultColorMap.h"
#include <algorithm>
#include <boost/tuple/tuple.hpp>

#if _OPENMP
# include <omp.h>
#endif

struct ConnectPairsResult
{
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

	ConnectPairsResult() :
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
		readMismatches(0)
	{}

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
			<< "max_depth_reverse" << "\t"
			<< "path_mismatches" << "\t"
			<< "read_mismatches" << "\n";
		return out;
	}

	friend std::ostream& operator <<(std::ostream& out,
		const ConnectPairsResult& o)
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
			<< o.maxDepthVisitedReverse << "\t"
			<< o.pathMismatches << "\t"
			<< o.readMismatches << "\n";

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
	bool fixErrors;
	bool longSearch;
	bool maskBases;

	ConnectPairsParams() :
		minMergedSeqLen(0),
		maxMergedSeqLen(1000),
		maxPaths(NO_LIMIT),
		maxBranches(NO_LIMIT),
		maxPathMismatches(NO_LIMIT),
		maxReadMismatches(NO_LIMIT),
		fixErrors(false),
		longSearch(false),
		maskBases(false) {}

};

static inline ConnectPairsResult connectPairs(
	unsigned k,
	const FastaRecord& read1,
	const FastaRecord& read2,
	const DBGBloom& g,
	const ConnectPairsParams& params)
{
	ConnectPairsResult result;

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

	unsigned startKmerPos = getStartKmerPos(k, read1, g,
			false, params.longSearch);

	unsigned goalKmerPos = getStartKmerPos(k, read2, g,
			true, params.longSearch);

	const FastaRecord* pRead1 = &read1;
	const FastaRecord* pRead2 = &read2;
	FastaRecord correctedRead1;
	FastaRecord correctedRead2;
	size_t unused;

	if (startKmerPos == NO_MATCH && params.fixErrors) {
		correctedRead1 = read1;
		if (correctSingleBaseError(g, k, correctedRead1, unused)) {
			startKmerPos = getStartKmerPos(k, correctedRead1, g,
				false, params.longSearch);
			assert(startKmerPos != NO_MATCH);
			pRead1 = &correctedRead1;
		}
	}

	if (goalKmerPos == NO_MATCH && params.fixErrors) {
		correctedRead2 = read2;
		if (correctSingleBaseError(g, k, correctedRead2, unused)) {
			goalKmerPos = getStartKmerPos(k, correctedRead2, g,
				true, params.longSearch);
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

	ConstrainedBidiBFSVisitor<DBGBloom> visitor(g, startKmer, goalKmer,
			params.maxPaths, minPathLen, maxPathLen, params.maxBranches);
	bidirectionalBFS(g, startKmer, goalKmer, visitor);

	std::vector< Path<Kmer> > paths;
	result.readNamePrefix = pRead1->id.substr(0, pRead1->id.find_last_of("/"));
	result.pathResult = visitor.pathsToGoal(paths);
	result.numNodesVisited = visitor.getNumNodesVisited();
	result.maxActiveBranches = visitor.getMaxActiveBranches();
	result.maxDepthVisitedForward = visitor.getMaxDepthVisited(FORWARD);
	result.maxDepthVisitedReverse = visitor.getMaxDepthVisited(REVERSE);

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

#if 0
# pragma omp critical(cerr)
	std::cerr << result;
#endif

	return result;
}

#endif
