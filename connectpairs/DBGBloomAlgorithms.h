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
#include "Graph/Path.h"
#include <climits>
#include <algorithm> // for std::max

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

#endif
