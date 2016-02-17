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
#include "lib/bloomfilter-2dfba08d120d7659e8c75cf5c501b3b9040e98cb/BloomFilter.hpp"

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>

#if _OPENMP
# include <omp.h>
#endif

namespace BloomDBG {

	/**
	 * Type for a vertex in the de Bruijn graph.
	 */
	typedef std::pair<MaskedKmer, RollingHash> Vertex;

	/**
	 * Parameters controlling assembly.
	 */
	struct AssemblyParams
	{
		/** Bloom filter size (in bits) */
		size_t bloomSize;

		/** minimum k-mer coverage threshold */
		unsigned minCov;

		/** path for output GraphViz file */
		string graphPath;

		/** approx genome size */
		size_t genomeSize;

		/** num Bloom filter hash functions */
		unsigned numHashes;

		/** the number of parallel threads. */
		unsigned threads;

		/** the size of a k-mer. */
		unsigned k;

		/** spaced seed */
		string spacedSeed;

		/** maximum length of branches to trim */
		unsigned trim;

		/** verbose level for progress messages */
		int verbose;

		/** Default constructor */
		AssemblyParams() : bloomSize(0), minCov(2), graphPath(), genomeSize(0),
			numHashes(1), threads(1), k(0), spacedSeed(),
			trim(std::numeric_limits<unsigned>::max()),
			verbose(0) {}

		/** Return true if all required members are initialized */
		bool initialized() const {
			return bloomSize > 0 && genomeSize > 0 && k > 0 &&
				!spacedSeed.empty() &&
				trim != std::numeric_limits<unsigned>::max();
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
	inline static void loadSeq(BF& bloom, const std::string& seq,
		const std::string& spacedSeed)
	{
		const unsigned k = bloom.getKmerSize();
		const unsigned numHashes = bloom.getHashNum();
		for (RollingHashIterator it(seq, k, numHashes, spacedSeed);
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
		const std::string& spacedSeed, bool verbose = false)
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
				loadSeq(bloom, buffer.at(j), spacedSeed);
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
	inline static bool allKmersInBloom(const Sequence& seq, const BloomT& bloom,
		const std::string& spacedSeed)
	{
		const unsigned k = bloom.getKmerSize();
		const unsigned numHashes = bloom.getHashNum();
		assert(seq.length() >= k);
		unsigned validKmers = 0;
		for (RollingHashIterator it(seq, k, numHashes, spacedSeed);
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
	inline static void addKmersToBloom(const Sequence& seq, BloomT& bloom,
		const std::string& spacedSeed)
	{
		const unsigned k = bloom.getKmerSize();
		const unsigned numHashes = bloom.getHashNum();
		for (RollingHashIterator it(seq, k, numHashes, spacedSeed);
			 it != RollingHashIterator::end(); ++it) {
			bloom.insert(*it);
		}
	}

	/**
	 * Translate a DNA sequence to an equivalent path in the
	 * de Bruijn graph.
	 */
	inline static Path<Vertex>
	seqToPath(const Sequence& seq, unsigned k, unsigned numHashes,
		const std::string spacedSeed)
	{
		Path<Vertex> path;
		assert(seq.length() >= k);
		for (RollingHashIterator it(seq, k, numHashes, spacedSeed);
			 it != RollingHashIterator::end(); ++it) {
			MaskedKmer kmer(it.kmer(), spacedSeed);
			path.push_back(Vertex(kmer, it.rollingHash()));
		}
		return path;
	}

	/**
	 * Translate a path in the de Bruijn graph to an equivalent
	 * DNA sequence.
	 */
	inline static Sequence pathToSeq(const Path<Vertex>& path, unsigned k,
		const std::string& spacedSeed)
	{
		assert(path.size() > 0);
		assert(k > 0);
		assert(spacedSeed.length() == k);

		Sequence seq;
		seq.resize(path.size() + k - 1, 'N');

		for (size_t i = 0; i < path.size(); ++i) {
			std::string kmer = path.at(i).first.str();
			for (size_t j = 0; j < k; ++j) {
				assert(spacedSeed.at(j) == '0' || spacedSeed.at(j) == '1');
				if (spacedSeed.at(j) == '1') {
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
	 * Extend a path left (REVERSE) or right (FORWARD) within the de Bruijn
	 * graph until either a branching point or a dead-end is encountered.
	 */
	template <typename GraphT>
	inline static bool extendPath(
		Path<typename boost::graph_traits<GraphT>::vertex_descriptor>& path,
		Direction dir, unsigned minBranchLen, const GraphT& graph)
	{
		unsigned origPathLen = path.size();
		assert(path.size() >= 1);

		/* Extend up to next branching point or dead end in DBG */
		extendPath(path, dir, graph, minBranchLen, NO_LIMIT);

		/* Return true if path was extended beyond orig length */
		return path.size() > origPathLen;
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
	 * Trim a sequence down to the longest contiguous subsequence
	 * of "good" k-mers.  If the sequence has length < k or contains
	 * no good k-mers, the trimmed sequence will be the empty string.
	 *
	 * @param seq the DNA sequence to be trimmed
	 * @param goodKmerSet Bloom filter containing "good" k-mers
	 * @param spacedSeed bitmap indicating positions to be ignored
	 * when hashing k-mers
	 */
	template <typename BloomT>
	static inline void trimSeq(Sequence& seq, const BloomT& goodKmerSet,
		const std::string& spacedSeed)
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
		for (RollingHashIterator it(seq, k, numHashes, spacedSeed);
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
	 * Trim first and/or last k-mer from a contig if they are
	 * branching k-mers that have already been included in other
	 * contigs.
	 */
	template <typename BloomT>
	void trimContig(Sequence& contig, const BloomT& assembledKmerSet,
		const std::string& spacedSeed)
	{
		const unsigned k = assembledKmerSet.getKmerSize();
		const unsigned numHashes = assembledKmerSet.getHashNum();
		std::vector<size_t> hashes;

		/* trim first k-mer */

		assert(contig.length() >= k);
		Sequence firstKmer = contig.substr(0, k);
		hashes = RollingHash(firstKmer, numHashes, k, spacedSeed).getHash();
		if (assembledKmerSet.contains(hashes)) {
			if (contig.length() == k) {
				contig.clear();
				return;
			}
			contig.erase(contig.begin());
		}

		/* trim last k-mer */

		assert(contig.length() >= k);
		Sequence lastKmer = contig.substr(contig.length() - k);
		hashes = RollingHash(lastKmer, numHashes, k, spacedSeed).getHash();
		if (assembledKmerSet.contains(hashes)) {
			if (contig.length() == k) {
				contig.clear();
				return;
			}
			contig.erase(contig.end() - 1);
		}
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
	 * @param spacedSeed bitmask indicating k-mer positions to ignore
	 * during hashing
	 * @param out output stream for contigs (FASTA)
	 * @param verbose set to true to print progress messages to
	 * STDERR
	 */
	template <typename BloomT>
	inline static void assemble(int argc, char** argv, const BloomT& goodKmerSet,
		const AssemblyParams& params, std::ostream& out)
	{
		assert(params.initialized());

		/* FASTA ID for next output contig */
		size_t contigID = 0;
		/* print progress message after processing this many reads */
		const unsigned progressStep = 1000;
		const unsigned k = goodKmerSet.getKmerSize();
		const unsigned numHashes = goodKmerSet.getHashNum();
		const std::string& spacedSeed = params.spacedSeed;
		/* k-mers in previously assembled contigs */
		BloomFilter assembledKmerSet(roundUpToMultiple(params.genomeSize,
			(size_t)64), numHashes, k);

		/* counters for progress messages */
		AssemblyCounters counters;

		/* Boost graph API over Bloom filter */
		RollingBloomDBG<BloomT> graph(goodKmerSet);

		unsigned minBranchLen = params.trim + 1;
		if (params.verbose)
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

			if (rec.seq.length() < k)
				skip = true;

			/* only extend error-free reads */
			if (!skip && !allKmersInBloom(rec.seq, goodKmerSet, spacedSeed))
				skip = true;

			/* skip reads in previously assembled regions */
			if (!skip && allKmersInBloom(rec.seq, assembledKmerSet, spacedSeed))
				skip = true;

			if (!skip) {

				std::cerr << "Extending read: " << rec.id << std::endl;

				/* convert sequence to DBG path */
				Path<Vertex> path = seqToPath(rec.seq, k, numHashes, spacedSeed);

				/* split path at branching points to prevent over-assembly */
				std::vector< Path<Vertex> > paths =
					splitPath(path, graph, minBranchLen);

				/*
				 * Extend first and last paths only, since
				 * internal path components are bounded by branching
				 * points.
				 */
				extendPath(paths.front(), REVERSE, minBranchLen, graph);
				extendPath(paths.back(), FORWARD, minBranchLen, graph);

				for(std::vector< Path<Vertex> >::iterator it = paths.begin();
					it != paths.end(); ++it) {
					/* convert DBG path back to sequence */
					Sequence seq = pathToSeq(*it, k, spacedSeed);
					/*
					 * check against assembledKmerSet again to prevent race
					 * condition. (Otherwise, the same contig may be
					 * generated multiple times.)
					 */
#pragma omp critical(out)
					if (!allKmersInBloom(seq, assembledKmerSet, spacedSeed)) {
						/*
						 * remove redundant branching k-mers at start/end
						 * of contig
						 */
						trimContig(seq, assembledKmerSet, spacedSeed);
						if (!seq.empty()) {
							assert(seq.length() >= k);
							addKmersToBloom(seq, assembledKmerSet, spacedSeed);
							FastaRecord contig;
							std::ostringstream id;
							id << contigID++;
							id << " read:" << rec.id;
							assert(id.good());
							contig.id = id.str();
							contig.seq = seq;
							out << contig;
							unsigned len = seq.length();
#pragma omp atomic
							counters.basesAssembled += len;
						}
					}
				}  /* for each split path */
#pragma omp atomic
				counters.readsExtended++;

			} /* if (!skip) */

#pragma omp atomic
			counters.readsProcessed++;
			if (params.verbose && counters.readsProcessed % progressStep == 0)
				printProgressMessage(counters);

		} /* for each read */

		assert(in.eof());
		if (opt::verbose) {
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
			m_out << '\t' << v.first.str() << ";\n";
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
			m_out << '\t' << u.first.str() << " -> "
				<< v.first.str() << ";\n";
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
	 * @param spacedSeed bitmask indicating which k-mer positions
	 * to ignore during hashing
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
		const std::string& spacedSeed = params.spacedSeed;

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
			trimSeq(seq, kmerSet, params.spacedSeed);
			if (seq.length() > 0) {

				/* BFS traversal in forward dir */
				Vertex start(MaskedKmer(seq.substr(0, k), spacedSeed),
					RollingHash(seq.substr(0, k), numHashes, k, spacedSeed));
				breadthFirstSearch(dbg, start, visitor, colorMap);

				/* BFS traversal in reverse dir */
				Sequence rcSeq = reverseComplement(seq);
				Vertex rcStart(MaskedKmer(rcSeq.substr(0, k), spacedSeed),
					RollingHash(rcSeq.substr(0, k), numHashes, k, spacedSeed));
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

} /* BloomDBG namespace */

#endif
