/*
 * MIBFQuerySupport.hpp
 *
 * Purpose: To provide support for complex classification
 *
 * Functions for most accurate classification and faster heuristic classification
 * are split into sections.
 *
 * Contains support objects intended to be private per thread (copied)
 *
 *
 *  Created on: Jun 6, 2018
 *      Author: justin
 */

#ifndef MIBFQUERYSUPPORT_HPP_
#define MIBFQUERYSUPPORT_HPP_

#include "MIBloomFilter.hpp"
//#include <set>
#include "vendor/ntHashIterator.hpp"
#include "vendor/stHashIterator.hpp"
#include <boost/math/distributions/binomial.hpp>

using namespace std;
using boost::math::binomial;

// T = T type, H = rolling hash itr
template<typename T>
class MIBFQuerySupport
{
  public:
	MIBFQuerySupport(
	    const MIBloomFilter<T>& miBF,
	    const vector<double>& perFrameProb,
	    unsigned extraCount,
	    unsigned extraFrameLimit,
	    unsigned maxMiss,
	    unsigned minCount,
	    bool bestHitAgree)
	  : m_miBF(miBF)
	  , m_perFrameProb(perFrameProb)
	  , m_extraCount(extraCount)
	  , m_extraFrameLimit(extraFrameLimit)
	  , m_maxMiss(maxMiss)
	  , m_minCount(minCount)
	  , m_bestHitAgree(bestHitAgree)
	  , m_satCount(0)
	  , m_evalCount(0)
	  , m_bestCounts({ 0, 0, 0, 0, 0, 0, 0 })
	  , m_secondBestNonSatFrameCount(0)
	  , m_rankPos(miBF.getHashNum())
	  , m_hits(miBF.getHashNum(), true)
	  , m_counts(vector<CountResult>(perFrameProb.size(), { 0, 0, 0, 0, 0, 0, 0 }))
	  , m_totalReads(0)
	{
		// this should always be a small array
		m_seenSet.reserve(miBF.getHashNum());
	}

	struct QueryResult
	{
		T id;
		uint16_t count;
		uint16_t nonSatCount;
		uint16_t totalCount;
		uint16_t totalNonSatCount;
		uint16_t nonSatFrameCount;
		uint16_t solidCount;
		double frameProb;
	};

	struct CountResult
	{
		uint16_t count;
		uint16_t nonSatCount;
		uint16_t totalCount;
		uint16_t totalNonSatCount;
		uint16_t nonSatFrameCount;
		uint16_t solidCount;
		size_t readCount; // determines if count should be reset
	};

	// For returning an empty result
	const vector<QueryResult>& emptyResult()
	{
		init();
		return m_signifResults;
	}

	/*
	 * totalTrials = number of possible trials that can be checked
	 */
	template<typename ITR>
	const vector<QueryResult>& query(ITR& itr, const vector<unsigned>& minCount)
	{
		init();

		unsigned extraFrame = 0;
		bool candidateFound = false;

		while (itr != itr.end() && !candidateFound) {
			candidateFound = updateCounts(itr, minCount, extraFrame);
			++itr;
		}
		summarizeCandiates();

		return m_signifResults;
	}

	template<typename ITR>
	const vector<QueryResult>& query(ITR& itr1, ITR& itr2, const vector<unsigned>& minCount)
	{
		init();

		unsigned extraFrame = 0;
		unsigned frameCount = 0;
		bool candidateFound = false;

		while ((itr1 != itr1.end() || itr2 != itr2.end()) && !candidateFound) {
			auto& itr =
			    frameCount % 2 == 0 && itr1 != itr1.end() ? itr1 : itr2 != itr2.end() ? itr2 : itr1;
			candidateFound = updateCounts(itr, minCount, extraFrame);
			++itr;
			++frameCount;
		}
		summarizeCandiates();

		return m_signifResults;
	}

	unsigned getSatCount() const { return m_satCount; }

	unsigned getEvalCount() const { return m_evalCount; }

	// debugging functions:

	void printAllCounts(const vector<string>& ids)
	{
		for (size_t i = 0; i < m_counts.size(); ++i) {
			if (m_counts[i].totalCount > 0) {
				cout << i << "\t" << ids[i] << "\t" << m_counts[i].nonSatFrameCount << "\t"
				     << m_counts[i].count << "\t" << m_counts[i].solidCount << "\t"
				     << m_counts[i].nonSatCount << "\t" << m_counts[i].totalNonSatCount << "\t"
				     << m_counts[i].totalCount << "\n";
			}
		}
	}

	/*
	 * Debugging
	 * Computes criteria used for judging a read consisting of:
	 * Position of matches
	 * Number of actually evaluated k-mers
	 * Return count of matching k-mers to set
	 */
	// TODO saturation not handle correctly
	inline vector<unsigned> getMatchSignature(
	    const string& seq,
	    unsigned& evaluatedSeeds,
	    vector<vector<pair<T, bool>>>& hitsPattern)
	{
		vector<unsigned> matchPos;
		matchPos.reserve(seq.size() - m_miBF.getKmerSize());

		if (m_miBF.getSeedValues().size() > 0) {
			stHashIterator itr(
			    seq, m_miBF.getSeedValues(), m_miBF.getHashNum(), m_miBF.getKmerSize());
			while (itr != itr.end()) {
				if (m_maxMiss >= m_miBF.atRank(*itr, m_rankPos, m_hits, m_maxMiss)) {
					vector<T> results = m_miBF.getData(m_rankPos);
					vector<pair<T, bool>> processedResults(results.size(), pair<T, bool>(0, false));
					for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
						if (m_hits[i]) {
							T tempResult = results[i];
							if (tempResult > MIBloomFilter<T>::s_mask) {
								processedResults[i] =
								    pair<T, bool>(tempResult & MIBloomFilter<T>::s_antiMask, true);
							} else {
								processedResults[i] =
								    pair<T, bool>(tempResult & MIBloomFilter<T>::s_antiMask, false);
							}
						}
					}
					matchPos.push_back(itr.pos());
					hitsPattern.push_back(processedResults);
				}
				++itr;
				++evaluatedSeeds;
			}
		} else {
			ntHashIterator itr(seq, m_miBF.getHashNum(), m_miBF.getKmerSize());
			while (itr != itr.end()) {
				if (m_miBF.atRank(*itr, m_rankPos)) {
					vector<T> results = m_miBF.getData(m_rankPos);
					vector<pair<T, bool>> processedResults(results.size(), pair<T, bool>(0, false));
					if (results.size() > 0) {
						for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
							T tempResult = results[i];
							if (tempResult > MIBloomFilter<T>::s_mask) {
								processedResults[i] =
								    pair<T, bool>(tempResult & MIBloomFilter<T>::s_antiMask, true);
							} else {
								processedResults[i] =
								    pair<T, bool>(tempResult & MIBloomFilter<T>::s_antiMask, false);
							}
						}
						matchPos.push_back(itr.pos());
						hitsPattern.push_back(processedResults);
					}
				}
				++itr;
				++evaluatedSeeds;
			}
		}
		return matchPos;
	}

  private:
	/*
	 * Sort in order of
	 * nonSatFrameCount
	 * count
	 * solidCount
	 * nonSatCount
	 * totalNonSatCount
	 * totalCount
	 * frameProb
	 */
	static inline bool sortCandidates(const QueryResult& a, const QueryResult& b)
	{
		return (
		    b.nonSatFrameCount == a.nonSatFrameCount
		        ? (b.count == a.count
		               ? (b.solidCount == a.solidCount
		                      ? (b.nonSatCount == a.nonSatCount
		                             ? (b.totalNonSatCount == a.totalNonSatCount
		                                    ? (b.totalCount == a.totalCount
		                                           ? (a.frameProb > b.frameProb)
		                                           : a.totalCount > b.totalCount)
		                                    : a.totalNonSatCount > b.totalNonSatCount)
		                             : a.nonSatCount > b.nonSatCount)
		                      : a.solidCount > b.solidCount)
		               : a.count > b.count)
		        : a.nonSatFrameCount > b.nonSatFrameCount);
	}

	//	static inline bool sortCandidates(const QueryResult &a,
	//			const QueryResult &b, unsigned extraCount ) {
	//		return (isRoughlyEqual(b.count, a.count, extraCount) ?
	//				(isRoughlyEqual(b.totalNonSatCount, a.totalNonSatCount, extraCount) ?
	//				(isRoughlyEqual(b.nonSatFrameCount, a.nonSatFrameCount, extraCount) ?
	//				(isRoughlyEqual(b.solidCount, a.solidCount, extraCount) ?
	//				(isRoughlyEqual(b.nonSatCount, a.nonSatCount, extraCount) ?
	//				(isRoughlyEqual(b.totalCount, a.totalCount, extraCount) ?
	//				(a.frameProb > b.frameProb) :
	//					a.totalCount > b.totalCount) :
	//					a.nonSatCount > b.nonSatCount) :
	//					a.solidCount > b.solidCount) :
	//					a.nonSatFrameCount > b.nonSatFrameCount) :
	//					a.totalNonSatCount > b.totalNonSatCount) :
	//					a.count > b.count);
	//	}

	//	static inline bool sortCandidates(const QueryResult &a,
	//			const QueryResult &b ) {
	//		return (compareStdErr(b.count, a.count) ?
	//				(compareStdErr(b.totalNonSatCount, a.totalNonSatCount) ?
	//				(compareStdErr(b.nonSatFrameCount, a.nonSatFrameCount) ?
	//				(compareStdErr(b.solidCount, a.solidCount) ?
	//				(compareStdErr(b.nonSatCount, a.nonSatCount) ?
	//				(compareStdErr(b.totalCount, a.totalCount) ?
	//				(a.frameProb > b.frameProb) :
	//					a.totalCount > b.totalCount) :
	//					a.nonSatCount > b.nonSatCount) :
	//					a.solidCount > b.solidCount) :
	//					a.nonSatFrameCount > b.nonSatFrameCount) :
	//					a.totalNonSatCount > b.totalNonSatCount) :
	//					a.count > b.count);
	//	}

	/*
	 * Returns true if considered roughly equal
	 */
	static inline bool isRoughlyEqual(unsigned a, unsigned b, unsigned extraCount)
	{
		if (a > b) {
			return a <= b + extraCount;
		}
		return b <= a + extraCount;
	}

	/*
	 * Returns true if considered roughly equal
	 */
	static inline bool compareStdErr(unsigned a, unsigned b)
	{
		double stderrA = sqrt(a);
		double stderrB = sqrt(b);
		if (a > b) {
			return (double(a) - stderrA) <= (double(b) + stderrB);
		}
		return (double(b) - stderrB) <= (double(a) + stderrA);
	}

	/*
	 * Returns true if considered roughly equal or b is larger
	 */
	inline bool compareStdErrLarger(unsigned a, unsigned b) const
	{
		double stderrA = sqrt(a) * m_extraCount;
		double stderrB = sqrt(b) * m_extraCount;
		return (double(a) - stderrA) <= (double(b) + stderrB);
	}

	/*
	 * Returns true if considered roughly equal
	 */
	bool isRoughlyEqual(const CountResult& a, const CountResult& b, unsigned extraCount) const
	{
		return (
		    isRoughlyEqual(b.count, a.count, extraCount) &&
		    isRoughlyEqual(b.totalNonSatCount, a.totalNonSatCount, extraCount) &&
		    isRoughlyEqual(b.nonSatFrameCount, a.nonSatFrameCount, extraCount) &&
		    isRoughlyEqual(b.solidCount, a.solidCount, extraCount) &&
		    isRoughlyEqual(b.nonSatCount, a.nonSatCount, extraCount) &&
		    isRoughlyEqual(b.totalCount, a.totalCount, extraCount));
	}

	/*
	 * Returns true if considered roughly equal
	 */
	bool isValid(const CountResult& a, const CountResult& b) const
	{
		return (
		    compareStdErr(b.count, a.count) ||
		    compareStdErr(b.totalNonSatCount, a.totalNonSatCount) ||
		    compareStdErr(b.nonSatFrameCount, a.nonSatFrameCount) ||
		    compareStdErr(b.solidCount, a.solidCount) ||
		    compareStdErr(b.nonSatCount, a.nonSatCount) ||
		    compareStdErr(b.totalCount, a.totalCount));
	}

	/*
	 * Returns true if considered roughly equal
	 */
	bool isRoughlyEqualOrLarger(const QueryResult& a, const QueryResult& b) const
	{
		return (
		    compareStdErrLarger(a.count, b.count) &&
		    compareStdErrLarger(a.totalNonSatCount, b.totalNonSatCount) &&
		    compareStdErrLarger(a.nonSatFrameCount, b.nonSatFrameCount) &&
		    compareStdErrLarger(a.solidCount, b.solidCount) &&
		    compareStdErrLarger(a.nonSatCount, b.nonSatCount) &&
		    compareStdErrLarger(a.totalCount, b.totalCount));
	}

	bool checkCountAgreement(QueryResult b, QueryResult a)
	{
		return (
		    b.nonSatFrameCount >= a.nonSatFrameCount && b.count >= a.count &&
		    b.solidCount >= a.solidCount && b.nonSatCount >= a.nonSatCount &&
		    b.totalNonSatCount >= a.totalNonSatCount && b.totalCount >= a.totalCount);
	}

	// contains reference to parent
	const MIBloomFilter<T>& m_miBF;
	const vector<double>& m_perFrameProb;

	// not references, but shared other objects or static variables
	const double m_extraCount;
	const unsigned m_extraFrameLimit;
	const unsigned m_maxMiss;
	const unsigned m_minCount;
	const bool m_bestHitAgree;
	//	const double m_rateSaturated;

	// resusable variables
	unsigned m_satCount;
	unsigned m_evalCount;

	// current bestCounts
	CountResult m_bestCounts;
	uint16_t m_secondBestNonSatFrameCount;

	// resusable objects
	vector<uint64_t> m_rankPos;
	vector<bool> m_hits;
	vector<QueryResult> m_signifResults;
	vector<CountResult> m_counts;
	vector<T> m_candidateMatches;
	vector<T> m_seenSet;

	// Number of reads processed by object
	size_t m_totalReads;

	bool
	updateCounts(const stHashIterator& itr, const vector<unsigned>& minCount, unsigned& extraFrame)
	{
		bool candidateFound = false;
		unsigned misses = m_miBF.atRank(*itr, m_rankPos, m_hits, m_maxMiss);
		if (misses <= m_maxMiss) {
			candidateFound = updatesCounts(minCount, extraFrame, misses);
		}
		return candidateFound;
	}

	bool
	updateCounts(const ntHashIterator& itr, const vector<unsigned>& minCount, unsigned& extraFrame)
	{
		bool candidateFound = false;
		if (m_miBF.atRank(*itr, m_rankPos)) {
			candidateFound = updatesCounts(minCount, extraFrame);
		}
		++m_evalCount;
		return candidateFound;
	}

	void init()
	{
		m_candidateMatches.clear();
		m_signifResults.clear();
		m_satCount = 0;
		m_evalCount = 0;
		m_bestCounts = { 0, 0, 0, 0, 0, 0, 0 };
		m_secondBestNonSatFrameCount = 0;
		++m_totalReads;
	}

	bool updatesCounts(const vector<unsigned>& minCount, unsigned& extraFrame, unsigned misses = 0)
	{
		m_seenSet.clear();
		unsigned satCount = 0;
		for (unsigned i = 0; i < m_miBF.getHashNum(); ++i) {
			if (m_hits[i]) {
				T resultRaw = m_miBF.getData(m_rankPos[i]);
				++m_evalCount;
				bool saturated = false;
				T result = resultRaw;

				// check for saturation
				if (result > m_miBF.s_mask) {
					result &= m_miBF.s_antiMask;
					saturated = true;
					satCount++;
					// detemines if count should be reset
					if (m_totalReads != m_counts[result].readCount) {
						m_counts[result] = { 0, 0, 0, 0, 0, 0, m_totalReads };
					}
				} else {
					if (m_totalReads != m_counts[result].readCount) {
						m_counts[result] = { 0, 0, 0, 0, 0, 0, m_totalReads };
					}
					++m_counts[result].totalNonSatCount;
				}
				++m_counts[result].totalCount;
				if (find(m_seenSet.begin(), m_seenSet.end(), resultRaw) == m_seenSet.end()) {
					// check for saturation
					if (saturated) {
						// if the non-saturated version has not been seen before
						if (find(m_seenSet.begin(), m_seenSet.end(), result) == m_seenSet.end()) {
							// check is count is exceeded
							++m_counts[result].count;
						}
					} else {
						++m_counts[result].nonSatCount;
						// check is count is exceeded
						++m_counts[result].count;
					}
					m_seenSet.push_back(resultRaw);
				}
			}
		}
		if (satCount == 0) {
			for (typename vector<T>::iterator itr = m_seenSet.begin(); itr != m_seenSet.end();
			     ++itr) {
				++m_counts[*itr].nonSatFrameCount;
				if (misses == 0) {
					++m_counts[*itr].solidCount;
				}
			}
		} else {
			++m_satCount;
		}
		for (typename vector<T>::iterator itr = m_seenSet.begin(); itr != m_seenSet.end(); ++itr) {
			T result = *itr;
			if (result > m_miBF.s_mask) {
				// if non-saturated version already exists
				if (find(m_seenSet.begin(), m_seenSet.end(), result & m_miBF.s_antiMask) !=
				    m_seenSet.end()) {
					continue;
				}
				result &= m_miBF.s_antiMask;
			}
			if (m_counts[result].count >= minCount[result]) {
				if (find(m_candidateMatches.begin(), m_candidateMatches.end(), result) ==
				    m_candidateMatches.end()) {
					m_candidateMatches.push_back(result);
				}
				updateMaxCounts(m_counts[result]);
			} else if (m_candidateMatches.size() && m_counts[result].count >= m_bestCounts.count) {
				if (find(m_candidateMatches.begin(), m_candidateMatches.end(), result) ==
				    m_candidateMatches.end()) {
					m_candidateMatches.push_back(result);
				}
				updateMaxCounts(m_counts[result]);
			}
		}
		if (compareStdErr(m_bestCounts.totalNonSatCount, m_secondBestNonSatFrameCount)) {
			extraFrame = 0;
		}
		if (m_bestCounts.nonSatFrameCount > m_secondBestNonSatFrameCount) {
			if (m_extraFrameLimit < extraFrame++) {
				return true;
			}
		}
		return false;
	}

	void updateMaxCounts(const CountResult& count)
	{
		if (count.nonSatFrameCount > m_bestCounts.nonSatFrameCount) {
			m_bestCounts.nonSatFrameCount = count.nonSatFrameCount;
		} else if (count.nonSatFrameCount > m_secondBestNonSatFrameCount) {
			m_secondBestNonSatFrameCount = count.nonSatFrameCount;
		}
		if (count.count > m_bestCounts.count) {
			m_bestCounts.count = count.count;
		}
		if (count.nonSatCount > m_bestCounts.nonSatCount) {
			m_bestCounts.nonSatCount = count.nonSatCount;
		}
		if (count.solidCount > m_bestCounts.solidCount) {
			m_bestCounts.solidCount = count.solidCount;
		}
		if (count.totalCount > m_bestCounts.totalCount) {
			m_bestCounts.totalCount = count.totalCount;
		}
		if (count.totalNonSatCount > m_bestCounts.totalNonSatCount) {
			m_bestCounts.totalNonSatCount = count.totalNonSatCount;
		}
	}

	double
	calcSat(unsigned evaluatedValues, double singleEventProbSaturted, unsigned saturatedCount)
	{
		double probSaturated = 0;
		if (saturatedCount) {
			binomial bin(evaluatedValues, singleEventProbSaturted);
			probSaturated = cdf(bin, saturatedCount - 1);
		}
		return probSaturated;
	}

	void summarizeCandiates()
	{
		if (m_candidateMatches.size() && m_minCount <= m_bestCounts.nonSatFrameCount) {
			vector<QueryResult> signifResults;
			for (typename vector<T>::const_iterator candidate = m_candidateMatches.begin();
			     candidate != m_candidateMatches.end();
			     candidate++) {
				const CountResult& resultCount = m_counts[*candidate];
				if (isValid(resultCount, m_bestCounts)) {
					QueryResult result;
					result.id = *candidate;
					result.count = resultCount.count;
					result.nonSatCount = resultCount.nonSatCount;
					result.totalCount = resultCount.totalCount;
					result.totalNonSatCount = resultCount.totalNonSatCount;
					result.nonSatFrameCount = resultCount.nonSatFrameCount;
					result.solidCount = resultCount.solidCount;
					result.frameProb = m_perFrameProb.at(*candidate);
					signifResults.push_back(result);
				}
			}
			if (signifResults.size() > 1) {
				sort(signifResults.begin(), signifResults.end(), sortCandidates);
				//				sort(signifResults.begin(), signifResults.end(),
				//						bind(sortCandidates, placeholders::_1, placeholders::_2,
				//								m_extraCount));
				for (typename vector<QueryResult>::iterator candidate = signifResults.begin();
				     candidate != signifResults.end();
				     candidate++) {
					if (isRoughlyEqualOrLarger(signifResults[0], *candidate)) {
						m_signifResults.push_back(*candidate);
					}
				}
				if (m_bestHitAgree && m_signifResults.size() >= 2 &&
				    !checkCountAgreement(m_signifResults[0], m_signifResults[1])) {
					m_signifResults.clear();
				}
			} else {
				m_signifResults.push_back(signifResults[0]);
			}
		}
	}
};

#endif /* MIBFQUERYSUPPORT_HPP_ */
