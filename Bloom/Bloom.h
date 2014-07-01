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
	};

	/** for verbose option to track loading progress */
	static const unsigned LOAD_PROGRESS_STEP = 100000;
	/** file format version number */
	static const unsigned BLOOM_VERSION = 2;
	/** I/O buffer size when reading/writing bloom filter files */
	static const unsigned long IO_BUFFER_SIZE = 32*1024;

	/**
	 * How to treat existing bits in the bloom filter when
	 * reading in new data.
	 */
	enum LoadType {
		LOAD_OVERWRITE,
		LOAD_UNION,
		LOAD_INTERSECT
	};

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
	inline static void loadFile(BF& bloomFilter, unsigned k, const std::string& path,
			bool verbose = false)
	{
		assert(!path.empty());
		if (verbose)
			std::cerr << "Reading `" << path << "'...\n";
		FastaReader in(path.c_str(), FastaReader::FOLD_CASE);
		uint64_t count = 0;
		const size_t CHUNK_SIZE = 1000;
#pragma omp parallel
		for (std::string seq[CHUNK_SIZE];;) {
			bool good = true;
			size_t i = 0;
#pragma omp critical(in)
			for (; good && i < CHUNK_SIZE; i++)
				good = in >> seq[i];
			if (!good)
				i--;
			if (i == 0)
				break;
			for (size_t j = 0; j < i; j++) {
#pragma omp atomic
				count++;
				loadSeq(bloomFilter, k, seq[j]);
			}
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
		for (size_t i = 0, j = 0; i < bytes;) {
			size_t writeSize = std::min(IO_BUFFER_SIZE, bytes - i);
			for (size_t k = 0; k < writeSize; k++) {
				buf[k] = 0;
				for (unsigned l = 0; l < 8; l++, j++) {
					buf[k] <<= 1;
					if (j < bits && bloomFilter[j]) {
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

		assert(in);
		assert(header.startBitPos < header.fullBloomSize);
		assert(header.endBitPos < header.fullBloomSize);
		assert(header.startBitPos <= header.endBitPos);

		return header;
	}

	/** Read the bloom filter bit array from a stream */
	template <typename BF>
	static void readData(BF& bloomFilter, const Bloom::FileHeader& header,
			std::istream& in, LoadType loadType = LOAD_OVERWRITE,
			unsigned shrinkFactor = 1)
	{

		// shrink factor allows building a smaller
		// bloom filter from a larger one

		size_t size = header.fullBloomSize;

		if (size % shrinkFactor != 0) {
			std::cerr << "error: the number of bits in the original bloom "
				"filter must be evenly divisible by the shrink factor (`"
				<< shrinkFactor << "')\n";
			exit(EXIT_FAILURE);
		}

		size /= shrinkFactor;

		if((loadType == LOAD_UNION || loadType == LOAD_INTERSECT)
			&& size != bloomFilter.size()) {
			std::cerr << "error: can't union/intersect two bloom filters "
				"with different sizes.\n";
			exit(EXIT_FAILURE);
		} else {
			bloomFilter.resize(size);
		}

		// read bit vector

		if (loadType == LOAD_OVERWRITE)
			bloomFilter.reset();

		size_t offset = header.startBitPos;
		size_t bits = header.endBitPos - header.startBitPos + 1;
		size_t bytes = (bits + 7) / 8;

		char buf[IO_BUFFER_SIZE];
		for (size_t i = 0, j = offset; i < bytes; ) {
			size_t readSize = std::min(IO_BUFFER_SIZE, bytes - i);
			in.read(buf, readSize);
			assert(in);
			for (size_t k = 0; k < readSize; k++) {
				for (unsigned l = 0; l < 8 && j < offset + bits; l++, j++) {
					bool bit = buf[k] & (1 << (7 - l));
					size_t index = j % size;
					switch (loadType)
					{
					case LOAD_OVERWRITE:
					case LOAD_UNION:
						bloomFilter[index] |= bit;
						break;
					case LOAD_INTERSECT:
						bloomFilter[index] &= bit;
						break;
					}
				}
			}
			i += readSize;
		}

	}

	/** Read a bloom filter from a stream */
	template <typename BF>
	static void read(BF& bloomFilter, std::istream& in,
			LoadType loadType = LOAD_OVERWRITE,
			unsigned shrinkFactor = 1)
	{
		FileHeader header = readHeader(in);
		readData(bloomFilter, header, in, loadType, shrinkFactor);
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
