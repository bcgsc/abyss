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

	/** for verbose option to track loading progress */
	static const unsigned LOAD_PROGRESS_STEP = 100000;
	/** file format version number */
	static const unsigned BLOOM_VERSION = 2;
	/** I/O buffer size when reading/writing bloom filter files */
	static const unsigned long IO_BUFFER_SIZE = 32*1024;


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

	/** Load a sequence file into a bloom filter */
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
			loadSeq(bloomFilter, k, seq);
		}
		assert(in.eof());
		if (verbose)
			std::cerr << "Loaded " << count << " reads into bloom filter\n";
	}

	/** Load a sequence (string) into a bloom filter */
	template <typename BF>
	inline static void loadSeq(BF& bloomFilter, unsigned k, const std::string& seq)
	{
		if (seq.size() < k)
			return;
		for (size_t i = 0; i < seq.size() - k + 1; ++i) {
			std::string kmer = seq.substr(i, k);
			size_t pos = kmer.find_last_not_of("ACGTacgt");
			if (pos == std::string::npos) {
				bloomFilter.insert(Kmer(kmer));
			} else
				i += pos;
		}
	}

	/** Write a bloom filter to a stream */
	template <typename BF>
	static void write(const BF& bloomFilter, size_t fullBloomSize,
		size_t startBitPos, size_t endBitPos, std::ostream& out)
	{

		// file header

		out << BLOOM_VERSION << '\n';
		assert(out);
		out << Kmer::length() << '\n';
		assert(out);
		out << fullBloomSize
			<< '\t' << startBitPos
			<< '\t' << endBitPos
			<< '\n';
		assert(out);

		// bloom filter bits

		size_t bits = endBitPos - startBitPos + 1;
		size_t bytes = (bits + 7) / 8;
		char buf[IO_BUFFER_SIZE];
		for (size_t i = 0, j = startBitPos; i < bytes;) {
			size_t writeSize = std::min(IO_BUFFER_SIZE, bytes - i);
			for (size_t k = 0; k < writeSize; k++) {
				buf[k] = 0;
				for (unsigned l = 0; l < 8; l++, j++) {
					buf[k] <<= 1;
					if (j <= endBitPos && bloomFilter[j]) {
						buf[k] |= 1;
					}
				}
			}
			out.write(buf, writeSize);
			assert(out);
			i += writeSize;
		}
	}

	/** Write a bloom filter to a stream */
	template <typename BF>
	static void write(const BF& bloomFilter, std::ostream& out)
	{
		Bloom::write(bloomFilter, bloomFilter.size(), 0,
			bloomFilter.size() - 1, out);
	}

	//TODO: Bloom filter calculation methods
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
