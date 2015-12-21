#ifndef BLOOM_DBG_H
#define BLOOM_DBG_H 1

#include "BloomDBG/RollingHashIterator.h"
#include "Common/Uncompress.h"
#include "DataLayer/FastaReader.h"
#include "Graph/Path.h"
#include "Graph/ExtendPath.h"
#include "Common/Kmer.h"
#include "BloomDBG/RollingHash.h"
#include "BloomDBG/RollingBloomDBG.h"
#include "Common/UnorderedSet.h"
#include "DataLayer/FastaConcat.h"
#include "lib/bloomfilter-9061f087d8714109b865415f2ac05e4796d0cd74/BloomFilter.hpp"

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

#if _OPENMP
# include <omp.h>
#endif

namespace BloomDBG {

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
						std::cerr << "Loaded " << readCount
							<< " reads into bloom filter\n";
				}
			}
		}
		assert(in.eof());
		if (verbose) {
			std::cerr << "Loaded " << readCount << " reads from `"
					  << path << "` into bloom filter\n";
		}
	}

	/**
	 * Return true if all of the k-mers in `seq` are contained in `bloom`
	 * and false otherwise.
	 */
	template <typename BloomT>
	inline static bool allKmersInBloom(const Sequence& seq,
		const BloomT& bloom, unsigned k, unsigned numHashes)
	{
		for (RollingHashIterator it(seq, k, numHashes);
			 it != RollingHashIterator::end(); ++it) {
			if (!bloom.contains(*it))
				return false;
		}
		return true;
	}

	/**
	 * Add all k-mers of a DNA sequence to a Bloom filter.
	 */
	template <typename BloomT>
	inline static void addKmersToBloom(const Sequence& seq,
		BloomT& bloom, unsigned k, unsigned numHashes)
	{
		for (RollingHashIterator it(seq, k, numHashes);
			 it != RollingHashIterator::end(); ++it) {
			bloom.insert(*it);
		}
	}

	/**
	 * Translate a DNA sequence to an equivalent path in the
	 * de Bruijn graph.
	 */
	inline static Path< std::pair<Kmer, RollingHash> >
	seqToPath(const Sequence& seq, unsigned k, unsigned numHashes)
	{
		typedef std::pair<Kmer, RollingHash> V;
		Path<V> path;
		assert(seq.length() >= k);
		std::string kmer0 = seq.substr(0, k);
		for (RollingHashIterator it(seq, k, numHashes);
			 it != RollingHashIterator::end(); ++it) {
			Kmer kmer(it.kmer());
			path.push_back(V(kmer, it.rollingHash()));
		}
		return path;
	}

	/**
	 * Translate a path in the de Bruijn graph to an equivalent
	 * DNA sequence.
	 */
	inline static Sequence pathToSeq(
		const Path< std::pair<Kmer, RollingHash> >& path)
	{
		typedef std::pair<Kmer, RollingHash> V;
		Sequence seq;
		seq.reserve(path.size());
		assert(path.size() > 0);
		Path<V>::const_iterator it = path.begin();
		assert(it != path.end());
		seq.append(it->first.str());
		for (++it; it != path.end(); ++it)
			seq.append(1, it->first.getLastBaseChar());
		return seq;
	}

	/**
	 * Extend a path left and right within the de Bruijn graph until
	 * either a branching point or a dead-end is encountered.
	 */
	template <typename GraphT>
	inline static void extendPath(
	Path<typename boost::graph_traits<GraphT>::vertex_descriptor>& path,
		const GraphT& graph, unsigned minBranchLen)
	{
		typedef std::pair<Kmer, RollingHash> V;
		unsigned chopLen;

		/* track visited vertices to detect cycles */
		unordered_set<V> visited;

		/*
		 * extend path right
		 *
		 * note: start extending a few k-mers before the end, to reduce
		 * the chance of hitting a dead-end at a Bloom filter false positive
		 */
		chopLen = std::min(path.size() - 1, (size_t)minBranchLen);
		path.erase(path.end() - chopLen, path.end());
		extendPath(path, FORWARD, graph, visited, minBranchLen, NO_LIMIT);

		/*
		 * extend path left
		 *
		 * note: start extending a few k-mers after the start, to reduce
		 * the chance of hitting a dead-end at a Bloom filter false positive
		 */
		chopLen = std::min(path.size() - 1, (size_t)minBranchLen);
		path.erase(path.begin(), path.begin() + chopLen);
		extendPath(path, REVERSE, graph, visited, minBranchLen, NO_LIMIT);
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

		AssemblyCounters() : readsExtended(0), readsProcessed(0),
			basesAssembled(0) {}
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
		if (currentPath.size() > 1)
			splitPaths.push_back(currentPath);

		return splitPaths;
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
	inline static void assemble(int argc, char** argv, size_t genomeSize,
		const BloomT& goodKmerSet, std::ostream& out, bool verbose=false)
	{
		/* FASTA ID for next output contig */
		size_t contigID = 0;
		/* print progress message after processing this many reads */
		const unsigned progressStep = 1000;
		const unsigned k = goodKmerSet.k();
		const unsigned numHashes = goodKmerSet.numHashes();
		/* k-mers in previously assembled contigs */
		BloomFilter assembledKmerSet(roundUpToMultiple(genomeSize, (size_t)64),
			numHashes, k);

		/* counters for progress messages */
		AssemblyCounters counters;

		/* Boost graph API over Bloom filter */
		RollingBloomDBG<BloomT> graph(goodKmerSet);
		/* vertex type for de Bruijn graph */
		typedef std::pair<Kmer, RollingHash> V;

		/*
		 * calculate min length threshold for a "true branch"
		 * (not due to Bloom filter false positives)
		 */
		const double falseBranchProbability = 0.0001;
		const unsigned minBranchLen =
			(unsigned)ceil(log(falseBranchProbability)/log(goodKmerSet.FPR()));

		if (verbose)
			std::cerr << "Treating branches less than " << minBranchLen
				<< " k-mers as Bloom filter false positives" << std::endl;

		FastaConcat in(argv, argv + argc, FastaReader::FOLD_CASE);
#pragma omp parallel
		for (FastaRecord rec;;) {
			bool good;
#pragma omp critical(in)
			good = in >> rec;
			if (!good)
				break;

			bool skip = false;

			/* only extend error-free reads */
			if (!allKmersInBloom(rec.seq, goodKmerSet, k, numHashes))
				skip = true;

			/* skip reads in previously assembled regions */
			if (allKmersInBloom(rec.seq, assembledKmerSet, k, numHashes))
				skip = true;

			if (!skip) {

				/* convert sequence to DBG path */
				Path<V> path = seqToPath(rec.seq, k, numHashes);

				/* split path at branching points to prevent over-assembly */
				std::vector< Path<V> > paths =
					splitPath(path, graph, minBranchLen);
				for(std::vector< Path<V> >::iterator it = paths.begin();
					it != paths.end(); ++it) {
					/* extend first and last paths only */
					if (it == paths.begin() || it == paths.end()-1)
						extendPath(*it, graph, minBranchLen);
					/* convert DBG path back to sequence */
					Sequence seq = pathToSeq(*it);
					/*
					 * check against assembledKmerSet again to prevent race
					 * condition. Otherwise, the same contig may be
					 * generated multiple times.
					 */
#pragma omp critical(out)
					if (!allKmersInBloom(seq, assembledKmerSet, k, numHashes)) {
						addKmersToBloom(seq, assembledKmerSet, k, numHashes);
						FastaRecord contig;
						std::ostringstream id;
						id << contigID++;
						assert(id.good());
						contig.id = id.str();
						contig.seq = seq;
						out << contig;
						unsigned len = seq.length();
#pragma omp atomic
						counters.basesAssembled += len;
					}
				}  /* for each split path */
#pragma omp atomic
				counters.readsExtended++;

			} /* if (!skip) */

#pragma omp atomic
			counters.readsProcessed++;
			if (verbose && counters.readsProcessed % progressStep == 0)
				printProgressMessage(counters);

		} /* for each read */

		assert(in.eof());
		if (opt::verbose) {
			printProgressMessage(counters);
			std::cerr << "Assembly complete" << std::endl;
		}
	}

} /* BloomDBG namespace */

#endif
