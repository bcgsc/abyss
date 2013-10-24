/**
 * Algorithms for a de Bruijn Graph using a Bloom filter
 * Copyright 2013 Canada's Michael Smith Genome Science Centre
 */
#ifndef DBGBLOOMALGORITHMS_H
#define DBGBLOOMALGORITHMS_H 1

#include "Common/Kmer.h"
#include "DBGBloom.h"
#include "Common/StringUtil.h"
#include "DataLayer/FastaReader.h"
#include "Graph/DefaultColorMap.h"
#include "Graph/Path.h"
#include "Graph/BreadthFirstSearch.h"
#include "Graph/ConstrainedBFSVisitor.h"

#define SUPPRESS_UNUSED_WARNING(a) (void)a

static inline Sequence pathToSeq(Path<Kmer> path)
{
	Sequence seq;
	assert(path.size() > 0);
	seq.append(path[0].str());
	for (unsigned i = 1; i < path.size(); i++)
		seq.append(1, path[i].getLastBaseChar());
	return seq;
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
	SUPPRESS_UNUSED_WARNING(connectPairs);

	unsigned k = g.m_k;

	assert(isReadNamePair(read1.id, read2.id));

	if (read1.seq.length() < k || read2.seq.length() < k)
		return NO_PATH;

	std::string kmer1Str = read1.seq.substr(0, k);
	std::string kmer2Str = read2.seq.substr(0, k);

	// TODO: advance to next kmers in the reads instead of giving up

	if (kmer1Str.find_first_not_of("AGCTagct") != std::string::npos) {
		std::cerr << "failed to connect read pair: non-AGCT char in first kmer "
			<< "(read id = " << read1.id << ", kmer = " << kmer1Str << ")\n";
		return NO_PATH;
	}

	if (kmer2Str.find_first_not_of("AGCTagct") != std::string::npos) {
		std::cerr << "failed to connect read pair: non-AGCT char in first kmer "
			<< "(read id = " << read2.id << ", kmer = " << kmer2Str << ")\n";
		return NO_PATH;
	}

	// TODO: add option for mate pair orientation (RF)

	Kmer kmer1(kmer1Str);
	Kmer kmer2(kmer2Str);
	kmer2.reverseComplement();

	// TODO: advance to next kmers in the reads instead of giving up

	if (!graph_traits<DBGBloom>::vertex_exists(kmer1, g)) {
		std::cerr << "failed to connect read pair: bloom filter miss on first kmer "
			<< "(read id = " << read1.id << ", kmer = " << kmer1 << ")\n";
		return NO_PATH;
	}

	if (!graph_traits<DBGBloom>::vertex_exists(kmer2, g)) {
		std::cerr << "failed to connect read pair: bloom filter miss on last kmer "
			<< "(read id = " << read2.id << ", rc(kmer) = " << kmer2 << ")\n";
		return NO_PATH;
	}

	unsigned maxPathLen;
	if (maxMergedSeqLen == NO_LIMIT) {
		maxPathLen = NO_LIMIT;
	} else {
		assert(maxMergedSeqLen > 0);
		maxPathLen = maxMergedSeqLen - k + 1;
	}

	DefaultColorMap<DBGBloom> colorMap;
	// note: maxDepth param is maxPathLen - 1 because the start node is at depth 0 
	ConstrainedBFSVisitor<DBGBloom> visitor(kmer1, kmer2, 0, maxPathLen - 1, maxBranches, colorMap);
	breadthFirstSearch(g, kmer1, visitor, colorMap);

	std::vector< Path<Kmer> > pathsFound;
	PathSearchResult result = visitor.pathsToGoal(pathsFound, maxPaths);

	if (result == FOUND_PATH) {
		std::string mergedId = read1.id.substr(0, read1.id.find_last_of("/"));
		for (unsigned i = 0; i < pathsFound.size(); i++) {
			FastaRecord mergedSeq;
			mergedSeq.id = mergedId;
			mergedSeq.seq = pathToSeq(pathsFound[i]);
			mergedSeqs.push_back(mergedSeq);
		}
	}

	return result;
}

#endif
