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
#include "Graph/BreadthFirstSearch.h"
#include "Graph/ConstrainedBFSVisitor.h"
#include <climits>

#if _OPENMP
# include <omp.h>
#endif

#define NO_MATCH UINT_MAX

static inline Sequence pathToSeq(Path<Kmer> path)
{
	Sequence seq;
	assert(path.size() > 0);
	seq.append(path[0].str());
	for (unsigned i = 1; i < path.size(); i++)
		seq.append(1, path[i].getLastBaseChar());
	return seq;
}

static inline unsigned getStartKmerPos(const FastaRecord& read, const DBGBloom& g, bool rc = false)
{
	unsigned k = g.m_k;

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
				matchLength = 0;
				matchPosSet = false;
			}
		}
	}
	assert(maxMatchLength > 0);

	// return the kmer closest to the gap between read pairs
	return maxMatchPos + maxMatchLength - 1;
}

static inline PathSearchResult connectPairs(
	const FastaRecord& read1,
	const FastaRecord& read2,
	const DBGBloom& g,
	std::vector<FastaRecord>& mergedSeqs,
	unsigned maxPaths = 2,
	unsigned maxMergedSeqLen = NO_LIMIT,
	unsigned maxBranches = NO_LIMIT)
{
	unsigned k = g.m_k;

	assert(isReadNamePair(read1.id, read2.id));

	if (read1.seq.length() < k || read2.seq.length() < k) {
#pragma omp critical(cerr)
		std::cerr 
			<< "failed to connect read pair: first or second read length is less than k"
			<< "(read1 = " << read1.id << ", read2 = " << read2.id  << ")\n";
		return NO_PATH;
	}

	unsigned kmer1Pos = getStartKmerPos(read1, g, false);
	unsigned kmer2Pos = getStartKmerPos(read2, g, true);

	if (kmer1Pos == NO_MATCH || kmer2Pos == NO_MATCH) {
#pragma omp critical(cerr)
		std::cerr 
			<< "failed to connect read pair: couldn't find bloom filter match in first or second read "
			<< "(read1 = " << read1.id << ", read2 = " << read2.id  << ")\n";
		return NO_PATH;
	}

	Kmer kmer1(read1.seq.substr(kmer1Pos, k));
	Kmer kmer2(read2.seq.substr(kmer2Pos, k));
	kmer2.reverseComplement();

	unsigned maxPathLen = NO_LIMIT;
	if (maxMergedSeqLen != NO_LIMIT) {
		maxPathLen = maxMergedSeqLen - k + 1 - kmer1Pos - kmer2Pos;
		// check for overflow
		assert(maxPathLen < maxMergedSeqLen);
	}

	DefaultColorMap<DBGBloom> colorMap;
	// note: maxDepth param is maxPathLen - 1 because the start node is at depth 0 
	ConstrainedBFSVisitor<DBGBloom> visitor(kmer1, kmer2, 0, maxPathLen - 1, maxBranches, colorMap);
	breadthFirstSearch(g, kmer1, visitor, colorMap);

	std::vector< Path<Kmer> > pathsFound;
	PathSearchResult result = visitor.pathsToGoal(pathsFound, maxPaths);

	if (result == FOUND_PATH) {
		std::string mergedId = read1.id.substr(0, read1.id.find_last_of("/"));
		std::string seqPrefix = read1.seq.substr(0, kmer1Pos);
		std::string seqSuffix = reverseComplement(read2.seq.substr(0, kmer2Pos));
		for (unsigned i = 0; i < pathsFound.size(); i++) {
			FastaRecord mergedSeq;
			mergedSeq.id = mergedId;
			mergedSeq.seq = seqPrefix + pathToSeq(pathsFound[i]) + seqSuffix;
			mergedSeqs.push_back(mergedSeq);
		}
	}

	return result;
}

#endif
