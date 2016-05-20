#ifndef BLOOM_DBG_H
#define BLOOM_DBG_H 1

#include "BloomDBG/RollingHashIterator.h"
#include "Common/Uncompress.h"
#include "Common/IOUtil.h"
#include "DataLayer/FastaReader.h"
#include "Graph/Path.h"
#include "Graph/ExtendPath.h"
#include "Graph/BreadthFirstSearch.h"
#include "BloomDBG/MaskedKmer.h"
#include "BloomDBG/RollingHash.h"
#include "BloomDBG/RollingBloomDBG.h"
#include "Common/UnorderedSet.h"
#include "DataLayer/FastaConcat.h"
#include "lib/bloomfilter/BloomFilter.hpp"

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <string>

#if _OPENMP
# include <omp.h>
#endif

namespace BloomDBG {

	/**
	 * Type for a vertex in the de Bruijn graph.
	 */
	typedef RollingBloomDBGVertex Vertex;

	/**
	 * Parameters controlling assembly.
	 */
	struct AssemblyParams
	{
		/** Bloom filter size (in bits) */
		size_t bloomSize;

		/** minimum k-mer coverage threshold */
		unsigned minCov;

		/** WIG track containing 0/1 for sufficient k-mer cov */
		std::string covTrackPath;

		/** path for output GraphViz file */
		string graphPath;

		/** num Bloom filter hash functions */
		unsigned numHashes;

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
		string spacedSeed;

		/** maximum length of branches to trim */
		unsigned trim;

		/** verbose level for progress messages */
		int verbose;

		/** output contigs path (empty string indicates STDOUT) */
		std::string outputPath;

		/** output path for trace file (-T) option */
		std::string tracePath;

		/** Default constructor */
		AssemblyParams() : bloomSize(0), minCov(2), graphPath(),
			numHashes(1), threads(1), k(0), K(0), qrSeedLen(0),
			spacedSeed(), trim(std::numeric_limits<unsigned>::max()),
			verbose(0), outputPath(), tracePath() {}

		/** Return true if all required members are initialized */
		bool initialized() const {
			return bloomSize > 0 && k > 0 &&
				trim != std::numeric_limits<unsigned>::max();
		}

		/** Reset all spaced seed params to their default values */
		void resetSpacedSeedParams() {
			spacedSeed.clear();
			K = 0;
			qrSeedLen = 0;
		}
	};

	/**
	 * Round up `num` to the nearest multiple of `base`.
	 */
	template <typename T>
	inline static T roundUpToMultiple(T num, T base)
	{
		if (base == 0)
			return num;
		T remainder = num % base;
		if (remainder == 0)
			return num;
		return num + base - remainder;
	}

	/**
	 * Load DNA sequence into Bloom filter using rolling hash.
	 *
	 * @param bloom target Bloom filter
	 * @param seq DNA sequence
	 */
	template <typename BF>
	inline static void loadSeq(BF& bloom, const std::string& seq)
	{
		const unsigned k = bloom.getKmerSize();
		const unsigned numHashes = bloom.getHashNum();
		for (RollingHashIterator it(seq, numHashes, k);
			it != RollingHashIterator::end(); ++it) {
			bloom.insert(*it);
		}
	}

	/**
	 * Load sequences contents of FASTA file into Bloom filter using
	 * rolling hash.
	 * @param bloom target Bloom filter
	 * @param path path to FASTA file
	 * @param verbose if true, print progress messages to STDERR
	 */
	template <typename BF>
	inline static void loadFile(BF& bloom, const std::string& path,
		bool verbose = false)
	{
		const size_t BUFFER_SIZE = 1000000;
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
				loadSeq(bloom, buffer.at(j));
				if (verbose)
#pragma omp critical(cerr)
				{
					readCount++;
					if (readCount % LOAD_PROGRESS_STEP == 0)
						std::cerr << "Loaded " << readCount
							<< " reads into Bloom filter\n";
				}
			}
		}
		assert(in.eof());
		if (verbose) {
			std::cerr << "Loaded " << readCount << " reads from `"
					  << path << "` into Bloom filter\n";
		}
	}

	/**
	 * Return true if all of the k-mers in `seq` are contained in `bloom`
	 * and false otherwise.
	 */
	template <typename BloomT>
	inline static bool allKmersInBloom(const Sequence& seq, const BloomT& bloom)
	{
		const unsigned k = bloom.getKmerSize();
		const unsigned numHashes = bloom.getHashNum();
		assert(seq.length() >= k);
		unsigned validKmers = 0;
		for (RollingHashIterator it(seq, numHashes, k);
			 it != RollingHashIterator::end(); ++it, ++validKmers) {
			if (!bloom.contains(*it))
				return false;
		}
		/* if we skipped over k-mers containing non-ACGT chars */
		if (validKmers < seq.length() - k + 1)
			return false;
		return true;
	}

	/**
	 * Add all k-mers of a DNA sequence to a Bloom filter.
	 */
	template <typename BloomT>
	inline static void addKmersToBloom(const Sequence& seq, BloomT& bloom)
	{
		const unsigned k = bloom.getKmerSize();
		const unsigned numHashes = bloom.getHashNum();
		for (RollingHashIterator it(seq, numHashes, k);
			 it != RollingHashIterator::end(); ++it) {
			bloom.insert(*it);
		}
	}

	/**
	 * Translate a DNA sequence to an equivalent path in the
	 * de Bruijn graph.
	 */
	inline static Path<Vertex>
	seqToPath(const Sequence& seq, unsigned k, unsigned numHashes)
	{
		Path<Vertex> path;
		assert(seq.length() >= k);
		for (RollingHashIterator it(seq, numHashes, k);
			 it != RollingHashIterator::end(); ++it) {
			path.push_back(Vertex(it.kmer().c_str(), it.rollingHash()));
		}
		return path;
	}

	/**
	 * Translate a path in the de Bruijn graph to an equivalent
	 * DNA sequence.
	 */
	inline static Sequence pathToSeq(const Path<Vertex>& path, unsigned k)
	{
		assert(path.size() > 0);
		assert(k > 0);

		const std::string& spacedSeed = MaskedKmer::mask();
		assert(spacedSeed.empty() || spacedSeed.length() == k);
		Sequence seq;
		seq.resize(path.size() + k - 1, 'N');

		for (size_t i = 0; i < path.size(); ++i) {
			std::string kmer(path.at(i).kmer().c_str());
			for (size_t j = 0; j < k; ++j) {
				if (spacedSeed.empty() || spacedSeed.at(j) == '1') {
					if (seq.at(i + j) != 'N' && seq.at(i + j) != kmer.at(j)) {
						std::cerr
							<< "warning: inconsistent DBG path detected "
							"at position " << i + j << ": "
							<< seq.substr(0, i + j)
							<< " (orig base: '" << seq.at(i + j) << "'"
							<< ", new base: '" << kmer.at(j) << "')"
							<< std::endl;
					}
					seq.at(i + j) = kmer.at(j);
				}
			}
		}

		return seq;
	}

	/**
	 * Results for the extension of a read segment.
	 * Each instance represents a row in the trace file generated
	 * by the '-T' option for abyss-bloom-dbg.
	 */
	struct SeqExtensionResult
	{
		/** FASTA ID for origin read */
		std::string readId;
		/**
		 * Index of this segment within the read. (Prior to extension,
		 * each read is split into segments at branching k-mers.)
		 */
		unsigned readSegmentId;
		/** Total number of segments after splitting the read */
		unsigned numReadSegments;
		/** True if leftwards sequence extension was attempted */
		bool extendedLeft;
		/** True if rightwards sequence extension was attempted */
		bool extendedRight;
		/** Result code for attempted left sequence extension (e.g. DEAD END) */
		PathExtensionResult leftExtensionResult;
		/** Result code for attempted left sequence extension (e.g. DEAD END) */
		PathExtensionResult rightExtensionResult;
		/** Original length of the read segment prior to extension */
		unsigned origLength;
		/** length of left extension (bp) */
		unsigned leftExtensionLength;
		/** length of right extension (bp) */
		unsigned rightExtensionLength;
		/** total length of extended sequence (bp) */
		unsigned extendedLength;
		/**
		 * True if the extended sequence was excluded from the output contigs
		 * because it was redundant. (An identical sequence was generated
		 * when extending a previous read.)
		 */
		bool redundantContig;
		/** Contig ID assigned to extended segment */
		size_t contigID;

		SeqExtensionResult() :
			readId(),
			readSegmentId(std::numeric_limits<unsigned>::max()),
			numReadSegments(std::numeric_limits<unsigned>::max()),
			extendedLeft(false),
			extendedRight(false),
			leftExtensionResult(DEAD_END),
			rightExtensionResult(DEAD_END),
			origLength(std::numeric_limits<unsigned>::max()),
			leftExtensionLength(std::numeric_limits<unsigned>::max()),
			rightExtensionLength(std::numeric_limits<unsigned>::max()),
			extendedLength(std::numeric_limits<unsigned>::max()),
			redundantContig(false),
			contigID(std::numeric_limits<size_t>::max()) {}

		bool initialized() const
		{
			return !readId.empty() &&
				readSegmentId != std::numeric_limits<unsigned>::max() &&
				numReadSegments != std::numeric_limits<unsigned>::max() &&
				origLength != std::numeric_limits<unsigned>::max() &&
				leftExtensionLength != std::numeric_limits<unsigned>::max() &&
				rightExtensionLength != std::numeric_limits<unsigned>::max() &&
				extendedLength != std::numeric_limits<unsigned>::max();
		}

		static std::ostream& printHeaders(std::ostream& out)
		{
			out << "read_id\t"
				<< "read_segment_id\t"
				<< "num_read_segments\t"
				<< "left_extension_result\t"
				<< "right_extension_result\t"
				<< "orig_length\t"
				<< "left_extension_len\t"
				<< "right_extension_len\t"
				<< "extended_length\t"
				<< "redundant_contig\t"
				<< "contig_id\n";
			return out;
		}

		friend std::ostream& operator <<(std::ostream& out,
			const SeqExtensionResult& o)
		{
			if (o.redundantContig) {
				out << o.readId << '\t'
					<< o.readSegmentId << '\t'
					<< o.numReadSegments << '\t'
					<< "-\t"
					<< "-\t"
					<< o.origLength << '\t'
					<< "-\t"
					<< "-\t"
					<< "-\t"
					<< "true" << '\t'
					<< "-\n";
			} else {
				out << o.readId << '\t'
					<< o.readSegmentId << '\t'
					<< o.numReadSegments << '\t';
				if (o.extendedLeft)
					out << pathExtensionResultStr(o.leftExtensionResult) << '\t';
				else
					out << "-\t";
				if (o.extendedRight)
					out << pathExtensionResultStr(o.rightExtensionResult) << '\t';
				else
					out << "-\t";
				out << o.origLength << '\t';
				if (o.extendedLeft)
					out << o.leftExtensionLength << '\t';
				else
					out << "-\t";
				if (o.extendedRight)
					out << o.rightExtensionLength << '\t';
				else
					out << "-\t";
				out << o.extendedLength << '\t'
					<< "false" << '\t'
					<< o.contigID << '\n';
			}
			return out;
		}
	};

	/**
	 * Extend a sequence left (REVERSE) or right (FORWARD) within the de Bruijn
	 * graph until either a branching point or a dead-end is encountered.
	 */
	template <typename GraphT>
	inline static PathExtensionResult extendSeq(Sequence& seq, Direction dir,
		unsigned k, unsigned numHashes, unsigned minBranchLen, const GraphT& graph)
	{
		assert(seq.length() >= k);

		/* Convert sequence to path in DBG */
		Path<Vertex> path = seqToPath(seq, k, numHashes);

		/* Extend path */
		ExtendPathParams params;
		params.trimLen = minBranchLen - 1;
		params.maxLen = NO_LIMIT;
		PathExtensionResult result =
			extendPath(path, dir, graph, params);

		/* Convert extended path back to sequence */
		Sequence extendedSeq = pathToSeq(path, k);

		/*
		 * If a spaced seed is in effect, short paths may result in
		 * sequences containing 'N's.  However, since we only extend
		 * "perfect reads", we can replace the 'N's with the correct
		 * bases by overlaying the seed sequence.
		 */
		if (dir == FORWARD) {
			overlaySeq(seq, extendedSeq, 0);
		} else {
			assert(dir == REVERSE);
			overlaySeq(seq, extendedSeq, extendedSeq.length() - seq.length());
		}

		/*
		 * Replace orig seq with extended version.
		 */
		seq = extendedSeq;

		/* Return true if sequence was successfully extended */
		return result;
	}


	/**
	 * Counters for tracking assembly statistics and producing
	 * progress messages.
	 */
	struct AssemblyCounters
	{
		size_t readsExtended;
		size_t readsProcessed;
		size_t basesAssembled;
		size_t contigID;

		AssemblyCounters() : readsExtended(0), readsProcessed(0),
			basesAssembled(0), contigID(0) {}
	};

	/** Print an intermediate progress message during assembly */
	void printProgressMessage(AssemblyCounters counters)
	{
#pragma omp critical(cerr)
		std::cerr
			<< "Extended " << counters.readsExtended
			<< " of " << counters.readsProcessed
			<< " reads (" << std::setprecision(3) << (float)100
			* counters.readsExtended / counters.readsProcessed
			<< "%), assembled " << counters.basesAssembled
			<< " bp so far" << std::endl;
	}

	/**
	 * Split a path at branching k-mers (degree > 2).
	 */
	template <typename GraphT>
	inline static std::vector<
		Path<typename boost::graph_traits<GraphT>::vertex_descriptor> >
	splitPath(const Path<typename boost::graph_traits<GraphT>::vertex_descriptor>& path,
		const GraphT& dbg, unsigned minBranchLen)
	{
		assert(path.size() > 0);

		typedef typename boost::graph_traits<GraphT>::vertex_descriptor V;
		typedef typename Path<V>::const_iterator PathIt;

		std::vector< Path<V> > splitPaths;
		Path<V> currentPath;
		for (PathIt it = path.begin(); it != path.end(); ++it) {
			currentPath.push_back(*it);
			unsigned inDegree =
				trueBranches(*it, REVERSE, dbg, minBranchLen).size();
			unsigned outDegree =
				trueBranches(*it, FORWARD, dbg, minBranchLen).size();
			if (inDegree > 1 || outDegree > 1) {
				/* we've hit a branching point -- end the current
				 * path and start a new one */
				splitPaths.push_back(currentPath);
				currentPath.clear();
				currentPath.push_back(*it);
			}
		}
		if (currentPath.size() > 1 || splitPaths.empty())
			splitPaths.push_back(currentPath);

		assert(splitPaths.size() >= 1);
		return splitPaths;
	}

	/**
	 * Split a sequence at branching k-mers (degree > 2).
	 * Branching k-mers are shared between the resulting sequence
	 * segments.
	 */
	template <typename GraphT>
	inline static std::vector<Sequence>
	splitSeq(const Sequence& seq, unsigned k, unsigned numHashes,
		const GraphT& dbg, unsigned minBranchLen)
	{
		assert(seq.length() >= k);

		typedef typename boost::graph_traits<GraphT>::vertex_descriptor V;
		typedef typename Path<V>::const_iterator PathIt;

		std::vector<Sequence> segments;
		Path<V> path = seqToPath(seq, k, numHashes);
		PathIt start = path.begin();
		PathIt end = path.begin();

		for (; end != path.end(); ++end) {
			std::vector<V> inBranches = trueBranches(*end, REVERSE, dbg,
				minBranchLen-1);
			unsigned inDegree = inBranches.size();
			/*
			 * Tricky: Include the read itself in the list of valid
			 * incoming branches, even if it is shorter than trimLen.
			 */
			if (end > path.begin() && std::find(inBranches.begin(),
				inBranches.end(), *(end - 1)) == inBranches.end()) {
				inDegree++;
			}
			std::vector<V> outBranches = trueBranches(*end, FORWARD, dbg,
				minBranchLen-1);
			unsigned outDegree = outBranches.size();
			/*
			 * Tricky: Include the read itself in the list of valid
			 * outgoing branches, even if it is shorter than trimLen.
			 */
			if (end < path.end() - 1 && std::find(outBranches.begin(),
				outBranches.end(), *(end + 1)) == outBranches.end()) {
				outDegree++;
			}
			if (inDegree > 1 || outDegree > 1) {
				/* we've hit a branching point -- end the current
				 * segment and start a new one */
				Sequence segment = seq.substr(start - path.begin(),
					end - start + k);
				segments.push_back(segment);
				start = end;
			}
		}
		if (segments.empty() || segments.back().length() > k) {
			Sequence segment = seq.substr(start - path.begin(),
				end - start + k);
			segments.push_back(segment);
		}

		assert(segments.size() >= 1);
		return segments;
	}

	/**
	 * Trim a sequence down to the longest contiguous subsequence
	 * of "good" k-mers.  If the sequence has length < k or contains
	 * no good k-mers, the trimmed sequence will be the empty string.
	 *
	 * @param seq the DNA sequence to be trimmed
	 * @param goodKmerSet Bloom filter containing "good" k-mers
	 */
	template <typename BloomT>
	static inline void trimSeq(Sequence& seq, const BloomT& goodKmerSet)
	{
		const unsigned k = goodKmerSet.getKmerSize();
		const unsigned numHashes = goodKmerSet.getHashNum();

		if (seq.length() < k) {
			seq.clear();
			return;
		}

		const unsigned UNSET = UINT_MAX;
		unsigned prevPos = UNSET;
		unsigned matchStart = UNSET;
		unsigned matchLen = 0;
		unsigned maxMatchStart = UNSET;
		unsigned maxMatchLen = 0;

		/* note: RollingHashIterator skips over k-mer
		 * positions with non-ACGT chars */
		for (RollingHashIterator it(seq, numHashes, k);
			it != RollingHashIterator::end(); prevPos=it.pos(),++it) {
			if (!goodKmerSet.contains(*it) ||
				(prevPos != UNSET && it.pos() - prevPos > 1)) {
				/* end any previous match */
				if (matchStart != UNSET && matchLen > maxMatchLen) {
					maxMatchLen = matchLen;
					maxMatchStart = matchStart;
				}
				matchStart = UNSET;
				matchLen = 0;
			}
			if (goodKmerSet.contains(*it)) {
				/* initiate or extend match */
				if (matchStart == UNSET)
					matchStart = it.pos();
				matchLen++;
			}
		}
		/* handles case last match extends to end of seq */
		if (matchStart != UNSET && matchLen > maxMatchLen) {
			maxMatchLen = matchLen;
			maxMatchStart = matchStart;
		}
		/* if there were no matching k-mers */
		if (maxMatchLen == 0) {
			seq.clear();
			return;
		}
		/* trim read down to longest matching subseq */
		seq = seq.substr(maxMatchStart, maxMatchLen + k - 1);
	}

	/**
	 * Ensure that branching k-mers are not repeated in the output
	 * contigs by selectively trimming contig ends.
	 *
	 * The idea is to keep a branch k-mer if the edge leading to it
	 * is unambigous. For example, in the diagram below the contig
	 * generated from the right side would include the branching k-mer
	 * k5, whereas the two contigs entering on the left would discard it:
	 *
	 * ...-k1-k2
	 *          \
	 *           k5-k6-...
	 *          /
	 * ...-k3-k4
	 *
	 * @param seq the contig to trim
	 * @param k k-mer size
	 * @param numHashes number of Bloom filter hash functions
	 * @param minBranchLen minimum length of a "true" branch (shorter
	 * branches are assumed to be caused by sequencing errors or
	 * Bloom filter false positives).
	 */
	template <typename GraphT>
		inline static void trimBranchKmers(Sequence& seq,
			unsigned k, unsigned numHashes, unsigned minBranchLen,
			const GraphT& dbg)
	{
		assert(seq.length() >= k);

		if (seq.length() == k)
			return;

		Sequence firstKmer = seq.substr(0, k);
		Vertex vFirst(firstKmer.c_str(), RollingHash(firstKmer, numHashes, k));
		unsigned outDegree = trueDegree(vFirst, FORWARD, dbg, minBranchLen - 1);
		if (outDegree > 1)
			seq.erase(0, 1);

		if (seq.length() == k)
			return;

		Sequence lastKmer = seq.substr(seq.length()-k);
		Vertex vLast(lastKmer.c_str(), RollingHash(lastKmer, numHashes, k));
		unsigned inDegree = trueDegree(vLast, REVERSE, dbg, minBranchLen - 1);
		if (inDegree > 1)
			seq.erase(seq.length()-1, 1);
	}

	/**
	 * Append a contig to the output FASTA stream.
	 */
    inline static void printContig(const Sequence& seq,
		size_t contigID, const std::string& readID, unsigned k,
		std::ostream& out)
	{
		assert(seq.length() >= k);

		FastaRecord contig;

		/* set FASTA id */
		std::ostringstream id;
		id << contigID;

		/* add FASTA comment indicating extended read id */
		std::ostringstream comment;
		comment << "read:" << readID;
		assert(id.good());
		contig.id = id.str();
		contig.comment = comment.str();

		/* set seq (in canonical orientation) */
		Sequence rcSeq = reverseComplement(seq);
		contig.seq = (seq < rcSeq) ? seq : rcSeq;

		/* output FASTQ record */
		out << contig;
		assert(out);
	}

	/**
	 * Trim contiguous stretches of previously-assembled k-mers from
	 * both ends of a contig.
	 *
	 * @param seq contig to be trimmed
	 * @param assembledKmerSet Bloom filter of k-mers from previously
	 * assembled contigs
	 */
	template <typename BloomT>
	inline static void trimContigOverlaps(Sequence &seq,
		const BloomT& assembledKmerSet)
	{
		const unsigned k = assembledKmerSet.getKmerSize();
		const unsigned numHashes = assembledKmerSet.getHashNum();

		/* trim previously assembled k-mers from start of sequence */
		RollingHashIterator fwd(seq, numHashes, k);
		for (; fwd != RollingHashIterator::end(); ++fwd) {
			if (!assembledKmerSet.contains(*fwd))
				break;
		}
		if (fwd.pos() > 0)
			seq.erase(0, fwd.pos());

		/* trim previously assembled k-mers from end of sequence */
		Sequence rcSeq = reverseComplement(seq);
		RollingHashIterator rev(rcSeq, numHashes, k);
		for (; rev != RollingHashIterator::end(); ++rev) {
			if (!assembledKmerSet.contains(*rev))
				break;
		}
		if (rev.pos() > 0)
			rcSeq.erase(0, rev.pos());

		/* flip seq back to original orientation */
		seq = reverseComplement(rcSeq);

		assert(seq.length() >= k);
	}

	/**
	 * Split a read at branching points in the de Bruijn graph and
	 * then extend each segment left and right, up to the next
	 * branching point or dead end.
	 *
	 * @param read read to be assembled
	 * @param dbg Boost graph interface to de Bruijn graph
	 * @param assembledKmerSet Bloom filter containing k-mers of
	 * previously assembled contigs
	 * @param params command line options for the assembly
	 * (e.g. k-mer coverage threshold)
	 * @param counters counter variables used for generating assembly
	 * progress messages.
	 * @param out output stream for contigs
	 * @param traceOut output stream for trace file (-T option)
	 */
	template <typename GraphT, typename BloomT>
	inline static void extendRead(const FastaRecord& read,
		const GraphT& dbg, BloomT& assembledKmerSet,
		const AssemblyParams& params, AssemblyCounters& counters,
		std::ostream& out, std::ostream& traceOut)
	{
		const unsigned k = params.k;
		const unsigned numHashes = params.numHashes;
		const unsigned minBranchLen = params.trim + 1;

		if (params.verbose >= 2) {
#pragma omp critical(cerr)
			std::cerr << "Extending read: " << read.id << std::endl;
		}

		/* split read at branching points (prevents over-assembly) */
		std::vector<Sequence> segments = splitSeq(read.seq, k,
			numHashes, dbg, minBranchLen);

		for (std::vector<Sequence>::iterator it = segments.begin();
			 it != segments.end(); ++it) {

			Sequence& seq = *it;

			/*
			 * track results of sequence extension attempt for
			 * trace file ('-T' option).
			 */
			SeqExtensionResult traceResult;
			traceResult.readId = read.id;
			traceResult.readSegmentId = it - segments.begin() + 1;
			traceResult.numReadSegments = segments.size();
			traceResult.origLength = seq.length();
			traceResult.leftExtensionLength = 0;
			traceResult.rightExtensionLength = 0;
			traceResult.redundantContig = true;

			/*
			 * extend first and last segments only, since
			 * internal segments are bounded by branching
			 * points.
			 */
			if (it == segments.begin()) {
				traceResult.extendedLeft = true;
				traceResult.leftExtensionResult = extendSeq(seq,
					REVERSE, k, numHashes, minBranchLen, dbg);
				traceResult.leftExtensionLength =
					seq.length() - traceResult.origLength;
			}
			if (it == segments.end() - 1) {
				unsigned origLength = seq.length();
				traceResult.extendedRight = true;
				traceResult.rightExtensionResult = extendSeq(seq,
					FORWARD, k, numHashes, minBranchLen, dbg);
				traceResult.rightExtensionLength =
					seq.length() - origLength;
			}
			traceResult.extendedLength = seq.length();

			/* ensure branching k-mers are included only once in output */
			trimBranchKmers(seq, k, numHashes, minBranchLen, dbg);

			/*
			 * check assembledKmerSet again to prevent race
			 * condition. (Otherwise, the same contig may be
			 * generated multiple times.)
			 */
#pragma omp critical(assembledKmerSet)
			if (!allKmersInBloom(seq, assembledKmerSet)) {

				/* trim previously assembled k-mers from both ends */
				trimContigOverlaps(seq, assembledKmerSet);

				/* mark remaining k-mers as assembled */
				addKmersToBloom(seq, assembledKmerSet);

				/* add contig to output FASTA */
				printContig(seq, counters.contigID, read.id, k, out);

				/* update counters / trace results */
				traceResult.redundantContig = false;
				traceResult.contigID = counters.contigID;
				counters.basesAssembled += seq.length();
				counters.contigID++;
			}

			/* trace file output ('-T' option) */
#pragma omp critical(traceOut)
			if (!params.tracePath.empty()) {
				assert(traceResult.initialized());
				traceOut << traceResult;
				assert_good(traceOut, params.tracePath);
			}

		}  /* for each read segment */
	}

	/**
	 * Perform a Bloom-filter-based de Bruijn graph assembly.
	 * Contigs are generated by extending reads left/right within
	 * the de Bruijn graph, up to the next branching point or dead end.
	 * Short branches due to Bloom filter false positives are
	 * ignored.
	 *
	 * @param argc number of input FASTA files
	 * @param argv array of input FASTA filenames
	 * @param genomeSize approx genome size
	 * @param goodKmerSet Bloom filter containing k-mers that
	 * occur more than once in the input data
	 * @param out output stream for contigs (FASTA)
	 * @param verbose set to true to print progress messages to
	 * STDERR
	 */
	template <typename BloomT>
	inline static void assemble(int argc, char** argv, const BloomT& goodKmerSet,
		const AssemblyParams& params, std::ostream& out)
	{
		assert(params.initialized());

		/* per-thread I/O buffer (size is in bases) */
		const size_t SEQ_BUFFER_SIZE = 1000000;

		/* print progress message after processing this many reads */
		const unsigned progressStep = 1000;
		const unsigned k = goodKmerSet.getKmerSize();

		/* trace file output ('-T' option) */
		std::ofstream traceOut;
		if (!params.tracePath.empty()) {
			traceOut.open(params.tracePath.c_str());
			assert_good(traceOut, params.tracePath);
			SeqExtensionResult::printHeaders(traceOut);
			assert_good(traceOut, params.tracePath);
		}

		/* k-mers in previously assembled contigs */
		BloomFilter assembledKmerSet(goodKmerSet.size(),
			goodKmerSet.getHashNum(), goodKmerSet.getKmerSize());
		/* counters for progress messages */
		AssemblyCounters counters;

		/* Boost graph API over Bloom filter */
		RollingBloomDBG<BloomT> graph(goodKmerSet);

		if (params.verbose)
			std::cerr << "Trimming branches " << params.trim
				<< " k-mers or shorter" << std::endl;

		FastaConcat in(argv, argv + argc, FastaReader::FOLD_CASE);
#pragma omp parallel
		for (std::vector<FastaRecord> buffer;;) {

			/* read sequences in batches to reduce I/O contention */
			buffer.clear();
			size_t bufferSize;
			bool good = true;
#pragma omp critical(in)
			for (bufferSize = 0; bufferSize < SEQ_BUFFER_SIZE;) {
				FastaRecord rec;
				good = in >> rec;
				if (!good)
					break;
				buffer.push_back(rec);
				bufferSize += rec.seq.length();
			}
			if (buffer.size() == 0)
				break;

			for (std::vector<FastaRecord>::iterator it = buffer.begin();
				 it != buffer.end(); ++it) {

				const FastaRecord& rec = *it;
				bool skip = false;

				/* we can't extend reads shorter than k */
				if (rec.seq.length() < k)
					skip = true;

				/* only extend error-free reads */
				if (!skip && !allKmersInBloom(rec.seq, goodKmerSet))
					skip = true;

				/* skip reads in previously assembled regions */
				if (!skip && allKmersInBloom(rec.seq, assembledKmerSet))
					skip = true;

				/* extend the read left and right within the DBG */
				if (!skip) {
					extendRead(rec, graph, assembledKmerSet, params,
						counters, out, traceOut);
#pragma omp atomic
					counters.readsExtended++;
				}

#pragma omp atomic
				counters.readsProcessed++;
				if (params.verbose && counters.readsProcessed % progressStep == 0)
					printProgressMessage(counters);

			} /* for each read */

		} /* for each batch of reads (parallel) */

		assert(in.eof());
		if (!params.tracePath.empty()) {
			traceOut.close();
			assert_good(traceOut, params.tracePath);
		}

		if (params.verbose) {
			printProgressMessage(counters);
			std::cerr << "Assembly complete" << std::endl;
		}
	}

	/**
	 * Visitor class that outputs visited nodes/edges in GraphViz format during
	 * a breadth first traversal. An instance of this class may be passed
	 * as an argument to the `breadthFirstSearch` function.
	 */
	template <typename GraphT>
	class GraphvizBFSVisitor
	{
		typedef typename boost::graph_traits<GraphT>::vertex_descriptor VertexT;
		typedef typename boost::graph_traits<GraphT>::edge_descriptor EdgeT;

	public:

		/** Constructor */
		GraphvizBFSVisitor(std::ostream& out) :
			m_out(out), m_nodesVisited(0), m_edgesVisited(0)
		{
			/* start directed graph (GraphViz) */
			m_out << "digraph g {\n";
		}

		/** Destructor */
		~GraphvizBFSVisitor()
		{
			/* end directed graph (GraphViz) */
			m_out << "}\n";
		}

		/** Invoked when a vertex is initialized */
		void initialize_vertex(const VertexT&, const GraphT&) {}

		/** Invoked when a vertex is visited for the first time */
		void discover_vertex(const VertexT& v, const GraphT&)
		{
			++m_nodesVisited;
			/* declare vertex (GraphViz) */
			m_out << '\t' << v.kmer().c_str() << ";\n";
		}

		/** Invoked each time a vertex is visited */
		void examine_vertex(const VertexT&, const GraphT&) {}

		/**
		 * Invoked when all of a vertex's outgoing edges have been
		 * traversed.
		 */
		void finish_vertex(const VertexT&, const GraphT&) {}

		/**
		 * Invoked when an edge is traversed. (Each edge
		 * in the graph is traversed exactly once.)
		 */
		void examine_edge(const EdgeT& e, const GraphT& g)
		{
			++m_edgesVisited;
			const VertexT& u = source(e, g);
			const VertexT& v = target(e, g);

			/* declare edge (GraphViz) */
			m_out << '\t' << u.kmer().c_str() << " -> "
				<< v.kmer().c_str() << ";\n";
		}

		/**
		 * Invoked when an edge is traversed to a "gray" vertex.
		 * A vertex is gray when some but not all of its outgoing edges
		 * have been traversed.
		 */
		void gray_target(const EdgeT&, const GraphT&) {}

		/**
		 * Invoked when an edge is traversed to a "black" vertex.
		 * A vertex is black when all of its outgoing edges have
		 * been traversed.
		 */
		void black_target(const EdgeT&, const GraphT&) {}

		/**
		 * Invoked when an edge is traversed to a "gray" or
		 * "black" vertex.
		 */
		void non_tree_edge(const EdgeT&, const GraphT&) {}

		/**
		 * Invoked when an edge is traversed to a "white" vertex.
		 * A vertex is a white if it is previously unvisited.
		 */
		void tree_edge(const EdgeT&, const GraphT&) {}

		/** Return number of distinct nodes visited */
		size_t getNumNodesVisited() const
		{
			return m_nodesVisited;
		}

		/** Get number of distinct edges visited */
		size_t getNumEdgesVisited() const
		{
			return m_edgesVisited;
		}

	protected:

		/** output stream for GraphViz serialization */
		std::ostream& m_out;
		/** number of nodes visited so far */
		size_t m_nodesVisited;
		/** number of edges visited so far */
		size_t m_edgesVisited;
	};

	/**
	 * Output a GraphViz serialization of the de Bruijn graph
	 * using FASTA files and a Bloom filter as input.
	 *
	 * @param argc number of input FASTA files
	 * @param argv array of input FASTA filenames
	 * @param kmerSet Bloom filter containing valid k-mers
	 * @param out output stream for GraphViz serialization
	 * @param verbose prints progress messages to STDERR if true
	 */
	template <typename BloomT>
	static inline void outputGraph(int argc, char** argv,
		const BloomT& kmerSet, const AssemblyParams& params,
		std::ostream& out)
	{
		assert(params.initialized());

		typedef RollingBloomDBG<BloomT> GraphT;

		/* interval for progress messages */
		const unsigned progressStep = 1000;
		const unsigned k = kmerSet.getKmerSize();
		const unsigned numHashes = kmerSet.getHashNum();

		/* counter for progress messages */
		size_t readsProcessed = 0;

		/* Boost graph API over rolling hash Bloom filter */
		GraphT dbg(kmerSet);

		/* Marks visited nodes in breadth-first traversal */
		DefaultColorMap<GraphT> colorMap;

		/* BFS Visitor -- generates GraphViz output as nodes
		 * and edges are traversed. */
		GraphvizBFSVisitor<GraphT> visitor(out);

		if (params.verbose)
			std::cerr << "Generating GraphViz output..." << std::endl;

		FastaConcat in(argv, argv + argc, FastaReader::FOLD_CASE);
		for (FastaRecord rec;;) {
			bool good;
			good = in >> rec;
			if (!good)
				break;
			Sequence& seq = rec.seq;

			/* Trim down to longest subsequence of "good" k-mers */
			trimSeq(seq, kmerSet);
			if (seq.length() > 0) {

				/* BFS traversal in forward dir */
				std::string startKmer = seq.substr(0, k);
				Vertex start(startKmer.c_str(),
					RollingHash(startKmer, numHashes, k));
				breadthFirstSearch(dbg, start, visitor, colorMap);

				/* BFS traversal in reverse dir */
				Sequence rcSeq = reverseComplement(seq);
				std::string rcStartKmer = rcSeq.substr(0, k);
				Vertex rcStart(rcStartKmer.c_str(),
					RollingHash(rcStartKmer, numHashes, k));
				breadthFirstSearch(dbg, rcStart, visitor, colorMap);

			}

			if (++readsProcessed % progressStep == 0 && params.verbose) {
				std::cerr << "processed " << readsProcessed
					<< " (k-mers visited: " << visitor.getNumNodesVisited()
					<< ", edges visited: " << visitor.getNumEdgesVisited()
					<< ")" << std::endl;
			}
		}
		assert(in.eof());
		if (params.verbose) {
			std::cerr << "processed " << readsProcessed
				<< " reads (k-mers visited: " << visitor.getNumNodesVisited()
				<< ", edges visited: " << visitor.getNumEdgesVisited()
				<< ")" << std::endl;
			std::cerr <<  "GraphViz generation complete" << std::endl;
		}
	}

	/**
	 * Write a single block of a 'variableStep' WIG file.
	 *
	 * @param chr chromosome name
	 * @param start start coordinate of block
	 * @param length length of block
	 * @param val value of block
	 * @param out output stream for WIG file
	 * @param outPath path for output WIG file
	 */
	static inline void outputWigBlock(const std::string& chr, size_t start,
		size_t length, unsigned val, ostream& out, const std::string& outPath)
	{
		assert(length > 0);
		out << "variableStep chrom=" << chr
			<< " span=" << length << "\n";
		out << start << ' ' << val << '\n';
		assert_good(out, outPath);
	}

	/**
	 * Write a WIG file for a reference genome, using the values 0 and 1
	 * to indicate whether or not a given k-mer had sufficient coverage
	 * in the reads to exceed the minimum coverage threshold.
	 *
	 * @param goodKmerSet Bloom filter of k-mers that exceed the
	 * minimum coverage threshold
	 * @param params encapsulates all command line options for the
	 * assembly, including the reference genome and the output path
	 * for the WIG file.
	 */
    template <class BloomT>
	static inline void writeCovTrack(const BloomT& goodKmerSet,
		const AssemblyParams& params)
	{
		assert(!params.covTrackPath.empty());
		assert(!params.refPath.empty());

		const unsigned k = goodKmerSet.getKmerSize();
		const unsigned numHashes = goodKmerSet.getHashNum();

		std::ofstream covTrack(params.covTrackPath.c_str());
		assert_good(covTrack, params.covTrackPath);

		if (params.verbose)
			std::cerr << "Writing 0/1 k-mer coverage track for `"
				<< params.refPath << "` to `"
				<< params.covTrackPath << "`" << std::endl;

		FastaReader ref(params.refPath.c_str(), FastaReader::FOLD_CASE);
		for (FastaRecord rec; ref >> rec;) {
			std::string chr = rec.id;
			bool firstVal = true;
			size_t blockStart = 1;
			size_t blockLength = 0;
			uint8_t blockVal = 0;
			for (RollingHashIterator it(rec.seq, numHashes, k);
				 it != RollingHashIterator::end(); ++it) {
				uint8_t val = goodKmerSet.contains(*it) ? 1 : 0;
				if (firstVal) {
					firstVal = false;
					/* WIG standard uses 1-based coords */
					blockStart = it.pos() + 1;
					blockLength = 1;
					blockVal = val;
				} else if (val != blockVal) {
					assert(firstVal == false);
					outputWigBlock(chr, blockStart, blockLength, blockVal,
						covTrack, params.covTrackPath);
					/* WIG standard uses 1-based coords */
					blockStart = it.pos() + 1;
					blockLength = 1;
					blockVal = val;
				} else {
					blockLength++;
				}
			}
			/* output last block */
			if (blockLength > 0) {
				outputWigBlock(chr, blockStart, blockLength, blockVal,
					covTrack, params.covTrackPath);
			}
		}
		assert(ref.eof());

		assert_good(covTrack, params.covTrackPath);
		covTrack.close();
	}

} /* BloomDBG namespace */

#endif
