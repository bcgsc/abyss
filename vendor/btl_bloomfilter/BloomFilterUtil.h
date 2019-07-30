#ifndef BLOOM_FILTER_UTIL
#define BLOOM_FILTER_UTIL 1

#include "KmerBloomFilter.hpp"
#include "vendor/ntHashIterator.hpp"

using namespace std;

void
insertSeq(BloomFilter& bloom, const string& seq, unsigned hashNum, unsigned kmerSize)
{
	ntHashIterator itr(seq, hashNum, kmerSize);
	while (itr != itr.end()) {
		bloom.insert(*itr);
		++itr;
	}
}

// functions for calculations regarding bloomfilter
// todo: Tweak calculations as they are approximations and may not be 100% optimal
// see http://en.wikipedia.org/wiki/Bloom_filter

// Private functions
/*
 * Calculate FPR based on hash functions, size and number of entries
 * see http://en.wikipedia.org/wiki/Bloom_filter
 */
double
calcApproxFPR(size_t size, size_t numEntr, unsigned hashFunctNum)
{
	return pow(
	    1.0 - pow(1.0 - 1.0 / double(size), double(numEntr) * hashFunctNum), double(hashFunctNum));
}

/*
 * Calculates redundancy FPR
 */
double
calcRedunancyFPR(size_t size, size_t numEntr, unsigned hashFunctNum)
{
	double total = log(calcApproxFPR(size, 1, hashFunctNum));
	for (size_t i = 2; i < numEntr; ++i) {
		total = log(exp(total) + calcApproxFPR(size, i, hashFunctNum));
	}
	return exp(total) / numEntr;
}

#endif
