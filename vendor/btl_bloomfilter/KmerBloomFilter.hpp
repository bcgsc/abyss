/*
 *
 * BloomFilter.hpp
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef KMERBLOOMFILTER_H_
#define KMERBLOOMFILTER_H_
#include "BloomFilter.hpp"
#include "vendor/nthash.hpp"

using namespace std;

class KmerBloomFilter : public BloomFilter
{
  public:
	/*
	 * Default constructor.
	 */
	KmerBloomFilter()
	  : BloomFilter()
	{}

	/* De novo filter constructor.
	 *
	 * preconditions:
	 * filterSize must be a multiple of 64
	 *
	 * kmerSize refers to the number of bases the kmer has
	 */
	KmerBloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize)
	  : BloomFilter(filterSize, hashNum, kmerSize)
	{}

	KmerBloomFilter(const string& filterFilePath)
	  : BloomFilter(filterFilePath)
	{}

	using BloomFilter::contains;
	using BloomFilter::insert;

	/*
	 * Single pass filtering, computes hash values on the fly
	 */
	bool contains(const char* kmer) const
	{
		uint64_t hVal = NTC64(kmer, m_kmerSize);
		size_t normalizedValue = hVal % m_size;
		unsigned char bit = bitMask[normalizedValue % bitsPerChar];
		if ((m_filter[normalizedValue / bitsPerChar] & bit) == 0)
			return false;
		for (unsigned i = 1; i < m_hashNum; i++) {
			normalizedValue = NTE64(hVal, m_kmerSize, i) % m_size;
			unsigned char bit = bitMask[normalizedValue % bitsPerChar];
			if ((m_filter[normalizedValue / bitsPerChar] & bit) == 0)
				return false;
		}
		return true;
	}

	void insert(const char* kmer)
	{
		uint64_t hVal = NTC64(kmer, m_kmerSize);
		size_t normalizedValue = hVal % m_size;
		__sync_fetch_and_or(
		    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
		for (unsigned i = 1; i < m_hashNum; i++) {
			size_t normalizedValue = NTE64(hVal, m_kmerSize, i) % m_size;
			__sync_fetch_and_or(
			    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
		}
	}
};

#endif /* KMERBLOOMFILTER_H_ */
