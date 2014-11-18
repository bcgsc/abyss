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
#include "Common/IOUtil.h"
#include "DataLayer/FastaReader.h"
#include <iostream>
#include <vector>

#if _OPENMP
# include <omp.h>
#endif

namespace Bloom {

	typedef Kmer key_type;

	/** Header section of serialized bloom filters. */
	struct FileHeader {
		unsigned bloomVersion;
		unsigned k;
		size_t fullBloomSize;
		size_t startBitPos;
		size_t endBitPos;
		size_t hashSeed;
	};

	/** Print a progress message after loading this many seqs */
	static const unsigned LOAD_PROGRESS_STEP = 100000;
	/** file format version number */
	static const unsigned BLOOM_VERSION = 4;

	/** Return the hash value of this object. */
	inline static size_t hash(const key_type& key)
	{
		if (key.isCanonical())
			return hashmem(&key, sizeof key);

		key_type copy(key);
		copy.reverseComplement();
		return hashmem(&copy, sizeof copy, 0);
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
	inline static void loadSeq(BF& bloomFilter, unsigned k, const std::string& seq);

	/** Load a sequence file into a bloom filter */
	template <typename BF>
	inline static void loadFile(BF& bloomFilter, unsigned k, const std::string& path,
			bool verbose = false, size_t taskIOBufferSize = 100000)
	{
		assert(!path.empty());
		if (verbose)
			std::cerr << "Reading `" << path << "'...\n";
		FastaReader in(path.c_str(), FastaReader::FOLD_CASE);
		uint64_t count = 0;
#pragma omp parallel
		for (std::vector<std::string> buffer(taskIOBufferSize);;) {
			buffer.clear();
			size_t bufferSize = 0;
			bool good = true;
#pragma omp critical(in)
			for (; good && bufferSize < taskIOBufferSize;) {
				std::string seq;
				good = in >> seq;
				if (good) {
					buffer.push_back(seq);
					bufferSize += seq.length();
				}
			}
			if (buffer.size() == 0)
				break;
			for (size_t j = 0; j < buffer.size(); j++) {
				loadSeq(bloomFilter, k, buffer.at(j));
				if (verbose)
#pragma omp critical(cerr)
				{
					count++;
					if (count % LOAD_PROGRESS_STEP == 0)
						std::cerr << "Loaded " << count << " reads into bloom filter\n";
				}
			}
		}
		assert(in.eof());
		if (verbose) {
			std::cerr << "Loaded " << count << " reads from `"
				<< path << "` into bloom filter\n";
		}
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

	inline static void writeHeader(std::ostream& out, const FileHeader& header)
	{
		(void)writeHeader;

		out << BLOOM_VERSION << '\n';
		out << Kmer::length() << '\n';
		out << header.fullBloomSize
			<< '\t' << header.startBitPos
			<< '\t' << header.endBitPos
			<< '\n';
		out << header.hashSeed << '\n';
		assert(out);
	}

	FileHeader readHeader(std::istream& in)
	{
		FileHeader header;

		// read bloom filter file format version

		in >> header.bloomVersion >> expect("\n");
		assert(in);
		if (header.bloomVersion != BLOOM_VERSION) {
			std::cerr << "error: bloom filter version (`"
				<< header.bloomVersion << "'), does not match version required "
				"by this program (`" << BLOOM_VERSION << "').\n";
			exit(EXIT_FAILURE);
		}

		// read bloom filter k value

		in >> header.k >> expect("\n");
		assert(in);
		if (header.k != Kmer::length()) {
			std::cerr << "error: this program must be run with the same kmer "
				"size as the bloom filter being loaded (k="
				<< header.k << ").\n";
			exit(EXIT_FAILURE);
		}

		// read bloom filter dimensions

		in >> header.fullBloomSize
		   >> expect("\t") >> header.startBitPos
		   >> expect("\t") >> header.endBitPos
		   >> expect("\n");

		in >> header.hashSeed >> expect("\n");

		assert(in);
		assert(header.startBitPos < header.fullBloomSize);
		assert(header.endBitPos < header.fullBloomSize);
		assert(header.startBitPos <= header.endBitPos);

		return header;
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
