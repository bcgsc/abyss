/**
 * Algorithms for a de Bruijn Graph using a Bloom filter
 * Copyright 2013 Canada's Michael Smith Genome Science Centre
 */
#ifndef DBGBLOOMALGORITHMS_H
#define DBGBLOOMALGORITHMS_H 1

#include "Common/Kmer.h"
#include "Common/KmerIterator.h"
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
	bool rc = false, bool longSearch = false)
{
	if (read.seq.size() < k)
		return NO_MATCH;

	// build a vector indicating whether each kmer is a match
	// Note: the vector intentionally has an extra false element at
	// the end, for the second loop below.

	const std::string& seq = read.seq;
	std::vector<bool> match(seq.length() - k + 2, false);
	bool foundMatch = false;
	for (unsigned i = 0; i < seq.length() - k + 1; i++) {
		std::string kmerStr = seq.substr(i,k);
		size_t pos = kmerStr.find_first_not_of("AGCTagct");
		if (pos != std::string::npos) {
			i += pos;
			continue;
		}
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
			if ((longSearch && matchLength > maxMatchLength) ||
				(!longSearch && matchLength >= maxMatchLength)) {
				maxMatchPos = matchPos;
				maxMatchLength = matchLength;
			}
			matchLength = 0;
			matchPosSet = false;
		}
	}
	assert(maxMatchLength > 0);

	if (longSearch)
		return maxMatchPos;
	else
		return maxMatchPos + maxMatchLength - 1;
}

struct BaseChangeScore
{

	size_t m_pos;
	char m_base;
	unsigned m_score;

public:

	BaseChangeScore() : m_pos(0), m_base('N'), m_score(0) { }

	BaseChangeScore(size_t pos, char base, unsigned score)
		: m_pos(pos), m_base(base), m_score(score) { }

};

static inline bool correctSingleBaseError(
		const DBGBloom& g,
		unsigned k,
		FastaRecord& read,
		size_t& correctedPos,
		bool rc = false)
{
	if (read.seq.length() < k)
		return false;

	SUPPRESS_UNUSED_WARNING(correctedPos);

	const std::string bases = "AGCT";
	const size_t minScore = 3;
	std::vector<BaseChangeScore> scores;

	for (size_t i = 0; i < read.seq.length(); i++) {

		size_t overlapStart = std::max((int)(i - k + 1), 0);
		size_t overlapEnd = std::min(i + k - 1, read.seq.length() - 1);
		assert(overlapStart < overlapEnd);
		Sequence overlapStr = read.seq.substr(overlapStart, overlapEnd - overlapStart + 1);
		size_t changePos = i - overlapStart;

		for (size_t j = 0; j < bases.size(); j++) {
			if (read.seq[i] == bases[j])
				continue;
			overlapStr[changePos] = bases[j];
			size_t score = 0;
			for (KmerIterator it(overlapStr, k, rc); it != KmerIterator::end(); it++) {
				if (graph_traits<DBGBloom>::vertex_exists(*it, g))
					score++;
			}
			if (score > minScore)
				scores.push_back(BaseChangeScore(i, bases[j], score));
		}

	}

	if (scores.size() == 0)
		return false;

	BaseChangeScore bestScore;
	bool bestScoreSet = false;
	for (size_t i = 0; i < scores.size(); i++) {
		if (!bestScoreSet || scores[i].m_score > bestScore.m_score) {
			bestScore = scores[i];
			bestScoreSet = true;
		}
	}

	correctedPos = bestScore.m_pos;
	read.seq[correctedPos] = bestScore.m_base;

	return true;
}


static inline SearchResult connectPairs(
	unsigned k,
	const FastaRecord& read1,
	const FastaRecord& read2,
	const DBGBloom& g,
	unsigned maxPaths = 2,
	unsigned minMergedSeqLen = 0,
	unsigned maxMergedSeqLen = NO_LIMIT,
	unsigned maxBranches = NO_LIMIT,
	bool fixErrors = false,
	bool longSearch = false)
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

	unsigned startKmerPos = getStartKmerPos(k, read1, g, false, longSearch);
	unsigned goalKmerPos = getStartKmerPos(k, read2, g, true, longSearch);

	const FastaRecord* pRead1 = &read1;
	const FastaRecord* pRead2 = &read2;
	FastaRecord correctedRead1;
	FastaRecord correctedRead2;
	size_t unused;

	if (startKmerPos == NO_MATCH && fixErrors) {
		correctedRead1 = read1;
		if (correctSingleBaseError(g, k, correctedRead1, unused)) {
			startKmerPos = getStartKmerPos(k, correctedRead1, g, false, longSearch);
			assert(startKmerPos != NO_MATCH);
			pRead1 = &correctedRead1;
		}
	}

	if (goalKmerPos == NO_MATCH && fixErrors) {
		correctedRead2 = read2;
		if (correctSingleBaseError(g, k, correctedRead2, unused, true)) {
			goalKmerPos = getStartKmerPos(k, correctedRead2, g, true, longSearch);
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

	unsigned maxPathLen = maxMergedSeqLen - k + 1 - startKmerPos - goalKmerPos;
	assert(maxPathLen <= maxMergedSeqLen - k + 1);

	unsigned minPathLen = (unsigned)std::max((int)0,
			(int)(minMergedSeqLen - k + 1 - startKmerPos - goalKmerPos));
	// do not allow merged seqs that are shorter than the reads
	minPathLen = std::max(minPathLen, (unsigned)std::max(
				pRead1->seq.length() - k + 1 - startKmerPos,
				pRead2->seq.length() - k + 1 - goalKmerPos));

	ConstrainedBidiBFSVisitor<DBGBloom> visitor(g, startKmer, goalKmer, maxPaths,
			minPathLen, maxPathLen, maxBranches);
	bidirectionalBFS(g, startKmer, goalKmer, visitor);

	std::vector< Path<Kmer> > paths;
	result.readNamePrefix = pRead1->id.substr(0, pRead1->id.find_last_of("/"));
	result.pathResult = visitor.pathsToGoal(paths);
	result.numNodesVisited = visitor.getNumNodesVisited();
	result.maxActiveBranches = visitor.getMaxActiveBranches();
	result.maxDepthVisitedForward = visitor.getMaxDepthVisited(FORWARD);
	result.maxDepthVisitedReverse = visitor.getMaxDepthVisited(REVERSE);

	if (result.pathResult == FOUND_PATH) {
		std::string seqPrefix = pRead1->seq.substr(0, startKmerPos);
		std::string seqSuffix = reverseComplement(pRead2->seq.substr(0, goalKmerPos));
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
