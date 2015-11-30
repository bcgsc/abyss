#ifndef BLOOM_DBG_H
#define BLOOM_DBG_H 1

#include "BloomDBG/RollingHashIterator.h"
#include "Common/Uncompress.h"
#include "DataLayer/FastaReader.h"
#include <string>
#include <iostream>

#if _OPENMP
# include <omp.h>
#endif

namespace BloomDBG {

	/**
	 * Load DNA sequence into Bloom filter using rolling hash.
	 *
	 * @param bloom target Bloom filter
	 * @param numHashes number of Bloom filter hash functions
	 * @param k k-mer size
	 * @param seq DNA sequence
	 */
	template <typename BF>
	inline static void loadSeq(BF& bloom, unsigned numHashes, unsigned k,
			const std::string& seq)
	{
		for (RollingHashIterator it(seq, k, numHashes);
			it != RollingHashIterator::end(); ++it) {
			bloom.insert(*it);
		}
	}

	/**
	 * Load sequences contents of FASTA file into Bloom filter using
	 * rolling hash.
	 * @param bloom target Bloom filter
	 * @param numHashes number of Bloom filter hash functions
	 * @param k k-mer size
	 * @param path path to FASTA file
	 * @param verbose if true, print progress messages to STDERR
	 */
	template <typename BF>
	inline static void loadFile(BF& bloom, unsigned numHashes, unsigned k,
		const std::string& path, bool verbose = false)
	{
		const size_t BUFFER_SIZE = 100000;
		const size_t LOAD_PROGRESS_STEP = 10000;

		assert(!path.empty());
		if (verbose)
			std::cerr << "Reading `" << path << "'..." << std::endl;

		FastaReader in(path.c_str(), FastaReader::FOLD_CASE);
		uint64_t readCount = 0;
#pragma omp parallel
		for (std::vector<std::string> buffer(BUFFER_SIZE);;) {
			buffer.clear();
			size_t bufferSize = 0;
			bool good = true;
#pragma omp critical(in)
			for (; good && bufferSize < BUFFER_SIZE;) {
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
				loadSeq(bloom, numHashes, k, buffer.at(j));
				if (verbose)
#pragma omp critical(cerr)
				{
					readCount++;
					if (readCount % LOAD_PROGRESS_STEP == 0)
						std::cerr << "Loaded " << readCount << " reads into bloom filter\n";
				}
			}
		}
		assert(in.eof());
		if (verbose) {
			std::cerr << "Loaded " << readCount << " reads from `"
					  << path << "` into bloom filter\n";
		}
	}

} /* BloomDBG namespace */

#endif
