#ifndef _ASSEMBLY_PARAMS_H_
#define _ASSEMBLY_PARAMS_H_

#include <string>
#include <iostream>
#include <limits>

namespace BloomDBG {

	/**
	 * Parameters controlling assembly.
	 */
	struct AssemblyParams
	{
		/** Bloom filter size (in bytes) */
		size_t bloomSize;

		/** Checkpoint frequency (reads processed per checkpoint) */
		size_t readsPerCheckpoint;

		/** Do not delete checkpoint files after a successful assembly */
		bool keepCheckpoint;

		/** Filename prefix for checkpoint files */
		std::string checkpointPathPrefix;

		/** minimum k-mer coverage threshold */
		unsigned minCov;

		/** WIG track containing 0/1 for sufficient k-mer cov */
		std::string covTrackPath;

		/** path for output GraphViz file */
		std::string graphPath;

		/** num Bloom filter hash functions */
		unsigned numHashes;

		/** input Bloom filter file (if empty, build Bloom filter from reads)*/
		std::string bloomPath;

		/** the number of parallel threads. */
		unsigned threads;

		/** the size of a k-mer. */
		unsigned k;

		/** the size of a single k-mer in a k-mer pair */
		unsigned K;

		/** reference genome */
		std::string refPath;

		/** Quadratic Residue (QR) seed length */
		unsigned qrSeedLen;

		/** spaced seed */
		std::string spacedSeed;

		/** maximum length of branches to trim */
		unsigned trim;

		/** verbose level for progress messages */
		int verbose;

		/** output contigs path (empty string indicates STDOUT) */
		std::string outputPath;

		/** output path for trace file (-T) option */
		std::string tracePath;

		/** Default constructor */
		AssemblyParams() : bloomSize(0),
			readsPerCheckpoint(std::numeric_limits<size_t>::max()),
			keepCheckpoint(false), checkpointPathPrefix("bloom-dbg-checkpoint"),
			minCov(2), graphPath(), numHashes(1), threads(1),
			k(0), K(0), qrSeedLen(0), spacedSeed(),
			trim(std::numeric_limits<unsigned>::max()),
			verbose(0), outputPath(), tracePath() {}

		/** Return true if all required members are initialized */
		bool initialized() const {
			return bloomSize > 0 && k > 0 &&
				trim != std::numeric_limits<unsigned>::max();
		}

		/** Return true if checkpoint creation is enabled */
		bool checkpointsEnabled() const {
			return readsPerCheckpoint != std::numeric_limits<size_t>::max();
		}

		/** Reset all spaced seed params to their default values */
		void resetSpacedSeedParams() {
			spacedSeed.clear();
			K = 0;
			qrSeedLen = 0;
		}

		/** Report current parameter values (for logging) */
		friend std::ostream& operator<<(std::ostream& out,
			const AssemblyParams& o)
		{
			out << "Assembly parameters:" << std::endl
				<< '\t' << "K-mer size (-k): " << o.k << std::endl
				<< '\t' << "K-mer coverage threshold (--kc): " << o.minCov << std::endl
				<< '\t' << "Max branch trim length (-t): " << o.trim << std::endl
				<< '\t' << "Bloom size in bytes (-b): " << o.bloomSize << std::endl
				<< '\t' << "Bloom hash functions (-H): " << o.numHashes << std::endl;

			if (o.K > 0)
				out << '\t' << "Spaced k-mer size (-K): " << o.K << std::endl;

			if (o.qrSeedLen > 0)
				out << '\t' << "Quadratic residue (QR) seed length (--qr-seed): "
					<< o.qrSeedLen << std::endl;

			return out;
		}
	};

} /* end of BloomDBG namespace */

#endif
