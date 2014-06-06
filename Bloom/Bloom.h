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

#include "Common/Kmer.h"
#include "Common/HashFunction.h"
#include "Common/Uncompress.h"
#include "DataLayer/FastaReader.h"
#include <iostream>

namespace Bloom {

	typedef Kmer key_type;

	//for verbose option to track loading progress
	static const unsigned LOAD_PROGRESS_STEP = 100000;

	/** Return the hash value of this object. */
	inline static size_t hash(const key_type& key)
	{
		if (key.isCanonical())
			return hashmem(&key, sizeof key);

		key_type copy(key);
		copy.reverseComplement();
		return hashmem(&copy, sizeof copy);
	}

	/** Return the hash value of this object given seed. */
	inline static size_t hash(const key_type& key, size_t seed)
	{
		if (key.isCanonical())
			return hashmem(&key, sizeof key, seed);

		key_type copy(key);
		copy.reverseComplement();
		return hashmem(&copy, sizeof copy, seed);
	}

	template <typename BF>
	inline static void loadFile(BF& bloomFilter, unsigned k, const std::string& path, bool verbose = false)
	{
		assert(!path.empty());
		if (verbose)
			std::cerr << "Reading `" << path << "'...\n";
		FastaReader in(path.c_str(), FastaReader::FOLD_CASE);
		uint64_t count = 0;
		for (std::string seq; in >> seq; count++) {
			if (verbose && count % LOAD_PROGRESS_STEP == 0)
				std::cerr << "Loaded " << count << " reads into bloom filter\n";
			bloomFilter.loadSeq(k, seq);
		}
		assert(in.eof());
		if (verbose)
			std::cerr << "Loaded " << count << " reads into bloom filter\n";
	}

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

};

#endif /* BLOOM_H_ */
