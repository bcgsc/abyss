#ifndef _CHECKPOINT_H_
#define _CHECKPOINT_H_

#include "BloomDBG/AssemblyCounters.h"
#include "BloomDBG/AssemblyParams.h"
#include "BloomDBG/AssemblyStreams.h"
#include "Common/IOUtil.h"
#include <cstdio>
#include <iostream>
#include <fstream>

namespace BloomDBG
{
	const std::string CHECKPOINT_FASTA_EXT = ".contigs.fa";
	const std::string CHECKPOINT_COUNTERS_EXT = ".counters.tsv";
	const std::string CHECKPOINT_BLOOM_DBG_EXT = ".dbg.bloom";
	const std::string CHECKPOINT_BLOOM_VISITED_EXT = ".visited.bloom";
	const std::string CHECKPOINT_TMP_EXT = ".tmp";

	/**
	 * Save the current assembly state to a set of checkpoint files,
	 * so that the assembly can be resumed from that point.
	 *
	 * @param dbg Bloom filter de Bruijn graph
	 * @param visitedKmerSet k-mers included in output contigs so far
	 * @param counters counters capturing assembly state
	 * (e.g. next input read index, next output contig index)
	 * @param params various assembly parameters (corresponding to
	 * command line opts)
	 */
	template <typename BloomDBGT, typename VisitedKmerSetT>
	static inline void createCheckpoint(const BloomDBGT& dbg,
		const VisitedKmerSetT& visitedKmerSet,
		const AssemblyCounters& counters,
		const AssemblyParams& params)
	{
		assert(params.checkpointsEnabled());
		assert(!params.checkpointPathPrefix.empty());
		std::string prefix = params.checkpointPathPrefix;

		if (params.verbose)
			std::cerr << "Writing checkpoint data..." << std::endl;

		/* write out Bloom filter de Bruijn graph */

		std::string dbgPath = prefix + CHECKPOINT_BLOOM_DBG_EXT;
		std::string dbgPathTmp = dbgPath + CHECKPOINT_TMP_EXT;

		if (params.verbose)
			std::cerr << '\t' << "Writing Bloom filter de Bruijn graph to `"
				<< dbgPathTmp << "'" << std::endl;

		std::ofstream dbgOut;
		dbgOut.open(dbgPathTmp.c_str());
		assert_good(dbgOut, dbgPathTmp);
		dbgOut << dbg;
		assert_good(dbgOut, dbgPathTmp);

		/* write out visited k-mers Bloom filter */

		std::string visitedPath = prefix + CHECKPOINT_BLOOM_VISITED_EXT;
		std::string visitedPathTmp = visitedPath + CHECKPOINT_TMP_EXT;

		if (params.verbose)
			std::cerr << '\t' << "Writing visited k-mers Bloom to `"
				<< visitedPathTmp << "'" << std::endl;

		std::ofstream visitedOut;
		visitedOut.open(visitedPathTmp.c_str());
		assert_good(visitedOut, visitedPathTmp);
		visitedOut << visitedKmerSet;
		assert_good(visitedOut, visitedPathTmp);

		/* write out index of next input read */

		std::string countersPath = prefix + CHECKPOINT_COUNTERS_EXT;
		std::string countersPathTmp = countersPath + CHECKPOINT_TMP_EXT;

		if (params.verbose)
			std::cerr << '\t' << "Writing assembly counters to `"
				<< countersPathTmp << "'" << std::endl;

		std::ofstream countersOut;
		countersOut.open(countersPathTmp.c_str());
		assert_good(countersOut, countersPathTmp);
		countersOut << counters;
		assert_good(countersOut, countersPathTmp);

		/* copy/move new checkpoint files on top of old ones */

		std::string fastaPath = prefix + CHECKPOINT_FASTA_EXT;
		std::string fastaPathTmp = fastaPath + CHECKPOINT_TMP_EXT;

		if (params.verbose)
			std::cerr << '\t' << "Copying `" << fastaPathTmp
				<< " to `" << fastaPath << "'" << std::endl;

		copyFile(fastaPathTmp, fastaPath);

		if (params.verbose)
			std::cerr << '\t' << "Moving `" << dbgPathTmp
				<< "' to `" << dbgPath << "'" << std::endl;

		if (rename(dbgPathTmp.c_str(), dbgPath.c_str()) != 0) {
			perror("Error renaming file");
			abort();
		}

		if (params.verbose)
			std::cerr << '\t' << "Moving `" << visitedPathTmp
				<< "' to `" << visitedPath << "'" << std::endl;

		if (rename(visitedPathTmp.c_str(), visitedPath.c_str()) != 0) {
			perror("Error renaming file");
			abort();
		}

		if (params.verbose)
			std::cerr << '\t' << "Moving `" << countersPathTmp
				<< "' to `" << countersPath << "'" << std::endl;

		if (rename(countersPathTmp.c_str(), countersPath.c_str()) != 0) {
			perror("Error renaming file");
			abort();
		}
	}

	/** Return true if checkpoint files exist and are readable */
	static inline bool checkpointExists(const AssemblyParams& params)
	{
		assert(params.checkpointsEnabled());
		assert(!params.checkpointPathPrefix.empty());
		std::string prefix = params.checkpointPathPrefix;

		std::string fastaPath = prefix + CHECKPOINT_FASTA_EXT;
		std::string dbgPath = prefix + CHECKPOINT_BLOOM_DBG_EXT;
		std::string visitedPath = prefix + CHECKPOINT_BLOOM_VISITED_EXT;
		std::string countersPath = prefix + CHECKPOINT_COUNTERS_EXT;

		return std::ifstream(fastaPath.c_str()).good()
			&& std::ifstream(dbgPath.c_str()).good()
			&& std::ifstream(visitedPath.c_str()).good()
			&& std::ifstream(countersPath.c_str()).good();
	}

	/**
	 * Restore assembly state from checkpoint files.
	 *
	 * @param dbg Bloom filter de Bruijn graph
	 * @param visitedKmerSet k-mers included in output contigs so far
	 * @param counters counters capturing assembly state
	 * (e.g. next input read index, next output contig index)
	 * @param params various assembly parameters (corresponding to
	 * @param out main output stream for assembled contigs
	 * command line opts)
	 */
	template <typename BloomDBGT, typename VisitedKmerSetT,
		typename InputReadStreamT>
	static inline void resumeFromCheckpoint(
		BloomDBGT& dbg, VisitedKmerSetT& visitedKmerSet,
		AssemblyCounters& counters, const AssemblyParams& params,
		AssemblyStreams<InputReadStreamT>& streams)
	{
		assert(checkpointExists(params));
		assert(params.checkpointsEnabled());
		assert(!params.checkpointPathPrefix.empty());
		std::string prefix = params.checkpointPathPrefix;

		if (params.verbose)
			std::cerr << "Resuming from last checkpoint..." << std::endl;

		/* load Bloom filter de Bruijn graph */

		std::string dbgPath = prefix + CHECKPOINT_BLOOM_DBG_EXT;
		if (params.verbose)
			std::cerr << '\t' << "Reading Bloom filter de Bruijn graph from `"
				<< dbgPath << "'" << std::endl;
		dbg.loadFilter(dbgPath);

		/* load visited k-mers Bloom filter */

		std::string visitedPath = prefix + CHECKPOINT_BLOOM_VISITED_EXT;
		if (params.verbose)
			std::cerr << '\t' << "Reading reading visited k-mers Bloom from `"
				<< visitedPath << "'" << std::endl;
		visitedKmerSet.loadFilter(visitedPath);

		/* load index for next input read */

		std::string countersPath = prefix + CHECKPOINT_COUNTERS_EXT;
		if (params.verbose)
			std::cerr << '\t' << "Reading index of next input read from `"
				<< countersPath << "'" << std::endl;
		std::ifstream countersIn(countersPath.c_str());
		assert_good(countersIn, countersPath);
		countersIn >> counters;
		assert_good(countersIn, countersPath);

		/* advance to previous position in input reads */

		if (params.verbose)
			std::cerr << '\t' << "Advancing to read index "
				<< counters.readsProcessed << " in input reads..." << std::endl;

		FastaRecord rec;
		for (size_t i = 0; i < counters.readsProcessed && streams.in; ++i)
			streams.in >> rec;
		assert (!streams.in.eof());

		/* restore previously assembled contigs */

		std::string fastaPath = prefix + CHECKPOINT_FASTA_EXT;
		std::string fastaPathTmp = fastaPath + CHECKPOINT_TMP_EXT;
		if (params.verbose)
			std::cerr << '\t' << "Copying `" << fastaPath
				<< "' to `" << fastaPathTmp << "'" << std::endl;
		copyFile(fastaPath, fastaPathTmp);

		if (params.verbose)
			std::cerr << '\t' << "Outputting previously assembled contigs "
				<< "from `" << fastaPath << "'" << std::endl;
		std::ifstream prevContigs(fastaPath.c_str());
		assert_good(prevContigs, fastaPath);
		streams.out << prevContigs.rdbuf();
		assert(streams.out);
	}

	/** Delete a file if it exists */
	static inline void removeFileIfExists(const std::string& path)
	{
		if (std::ifstream(path.c_str()).good()) {
			if (remove(path.c_str()) != 0) {
				perror("Error removing file");
				abort();
			}
		}
	}

	/** Remove all checkpoint-related files */
	static inline void removeCheckpointData(const AssemblyParams& params)
	{
		assert(params.checkpointsEnabled());
		assert(!params.checkpointPathPrefix.empty());
		std::string prefix = params.checkpointPathPrefix;

		if (params.verbose)
			std::cerr << "Removing checkpoint files..." << std::endl;

		/* remove Bloom filter de Bruijn graph file(s) */

		std::string dbgPath = prefix + CHECKPOINT_BLOOM_DBG_EXT;
		std::string dbgPathTmp = dbgPath + CHECKPOINT_TMP_EXT;
		removeFileIfExists(dbgPath);
		removeFileIfExists(dbgPathTmp);

		/* remove visited k-mers Bloom filter file(s) */

		std::string visitedPath = prefix + CHECKPOINT_BLOOM_VISITED_EXT;
		std::string visitedPathTmp = visitedPath + CHECKPOINT_TMP_EXT;
		removeFileIfExists(visitedPath);
		removeFileIfExists(visitedPathTmp);

		/* remove assembly counters file(s) */

		std::string countersPath = prefix + CHECKPOINT_COUNTERS_EXT;
		std::string countersPathTmp = countersPath + CHECKPOINT_TMP_EXT;
		removeFileIfExists(countersPath);
		removeFileIfExists(countersPathTmp);

		/* remove contigs FASTA file(s) */

		std::string fastaPath = prefix + CHECKPOINT_FASTA_EXT;
		std::string fastaPathTmp = fastaPath + CHECKPOINT_TMP_EXT;
		removeFileIfExists(fastaPath);
		removeFileIfExists(fastaPathTmp);
	}
}

#endif
