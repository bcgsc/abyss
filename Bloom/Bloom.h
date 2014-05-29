/*
 * Bloom.h
 *
 * This class was created when we removed the virtual class
 * BloomFilterBase.
 *
 * It contains a mix of components from BloomFilterBase as well as
 * additional functions and properties shared between Bloom Filter
 * classes.
 *
 * Contains useful algorithms for bloom filter calculations
 *
 *  Created on: May 28, 2014
 *      Author: cjustin
 */

#ifndef BLOOM_H_
#define BLOOM_H_

#include "Kmer.h"
#include "HashFunction.h"
#include "Uncompress.h"
#include "FastaReader.h"
#include <iostream>

class Bloom {
public:
	Bloom(size_t size, unsigned numHash);
	virtual ~Bloom();

	typedef Kmer key_type;

	/** Return the hash value of this object. */
	static size_t hash(const key_type& key)
	{
		if (key.isCanonical())
			return hashmem(&key, sizeof key);

		key_type copy(key);
		copy.reverseComplement();
		return hashmem(&copy, sizeof copy);
	}

	/** Return the hash value of this object given seed. */
	static size_t hash(const key_type& key, size_t seed)
	{
		if (key.isCanonical())
			return hashmem(&key, sizeof key, seed);

		key_type copy(key);
		copy.reverseComplement();
		return hashmem(&copy, sizeof copy, seed);
	}


private:
	const unsigned numHash;
	const size_t size; // Filter size as number of buckets
//	size_t uniqueEntries;
//	size_t redundantEntries;

};

//Bloom filter calculation methods
//static double calcApproxFPR(size_t numBucket, size_t numEntr,
//		unsigned numHash) const;
//static double calcRedunancyFPR(size_t numBucket, size_t numEntr,
//		unsigned numHash) const;
//static size_t calcOptimalSize(size_t numEle, double fpr) const;
//static size_t calcOptimalSize(size_t numEle, double fpr,
//		unsigned numHash) const;
//static unsigned calcOptimalNumHash(double fpr);
//static unsigned calcOptimalNumHash(size_t numBucket, size_t numEle);

#endif /* BLOOM_H_ */
