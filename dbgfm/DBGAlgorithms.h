/**
 * De Brujin Graph algorithms
 * Copyright 2013 Shaun Jackman
 */

#ifndef DBGALGORITHMS_H
#define DBGALGORITHMS_H 1

#include "Assembly/AssemblyAlgorithms.h"
#include "Common/Log.h"
#include "dbgfm/dbgfm.h"
#include <cassert>
#include <utility>
#include <limits.h>

namespace dbg {

using AssemblyAlgorithms::checkSeqContiguity;
using AssemblyAlgorithms::extendBranch;
using AssemblyAlgorithms::processLinearExtensionForBranch;

/** Generate the adjacency information for each sequence in the
 * collection.
 */
template<typename Graph>
void
generateAdjacency(Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::vertex_iterator Vit;
	typedef typename vertex_bundle_type<Graph>::type VP;

	Timer timer("GenerateAdjacency");

	size_t count = 0;
	size_t numBasesSet = 0;
	std::pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		VP& up = get(vertex_bundle, g, u);
		if (up.deleted())
			continue;

		if (++count % 1000000 == 0)
			logger(1) << "Finding adjacent k-mer: " << count << '\n';

		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			V v = u;
			v.shift(dir);
			for (unsigned i = 0; i < NUM_BASES; i++) {
				v.setLastBase(dir, i);
				if (vertex_exists(v, g)) {
					up.setBaseExtension(dir, i);
					numBasesSet++;
				}
			}
		}

		pumpNetwork(g);
	}

	if (numBasesSet > 0)
		logger(0) << "Added " << numBasesSet << " edges for "
			<< count << " vertices.\n";
}

/** Mark the specified vertex and its neighbours.
 * @return the number of marked edges
 */
template<typename Graph>
size_t
markNeighbours(Graph& g,
		typename graph_traits<Graph>::vertex_descriptor u,
		typename vertex_bundle_type<Graph>::type up,
		extDirection sense)
{
	size_t n = 0;
	Kmer v(u);
	v.shift(sense);
	for (unsigned i = 0; i < NUM_BASES; i++) {
		if (up.getExtension(sense).checkBase(i)) {
			++n;
			v.setLastBase(sense, i);
			if (sense == SENSE)
				put(vertex_mark_in, g, v, true);
			else
				put(vertex_mark_out, g, v, true);
		}
	}
	return n;
}

/** Mark ambiguous branches and branches from palindromes for removal.
 * @return the number of branches marked
 */
template<typename Graph>
size_t
markAmbiguous(Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::vertex_iterator Vit;
	typedef typename vertex_bundle_type<Graph>::type VP;

	Timer timer(__func__);
	size_t progress = 0;
	size_t countv = 0, counte = 0;

	std::pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		VP& up = get(vertex_bundle, g, u);
		if (up.deleted())
			continue;

		if (++progress % 1000000 == 0)
			logger(1) << "Splitting: " << progress << '\n';

		if (u.isPalindrome()) {
			countv += 2;
			up.mark();
			counte += markNeighbours(g, u, up, SENSE);
		} else {
			for (extDirection sense = SENSE;
					sense <= ANTISENSE; ++sense) {
				if (up.getExtension(sense).isAmbiguous()
						|| u.isPalindrome(sense)) {
					countv++;
					up.mark(sense);
					counte += markNeighbours(g, u, up, sense);
				}
			}
		}

		pumpNetwork(g);
	}
	logger(0) << "Marked " << counte << " edges of " << countv
		<< " ambiguous vertices." << std::endl;
	return countv;
}

/** Assemble a contig.
 * @return the number of k-mer below the coverage threshold
 */
template<typename Graph>
size_t
assembleContig(Graph& g, FastaWriter* writer,
		BranchRecord& branch, unsigned id)
{
	assert(!branch.isActive());
	assert(branch.getState() == BS_NOEXT
			|| branch.getState() == BS_AMBI_SAME
			|| branch.getState() == BS_AMBI_OPP);

	// Assemble the contig.
	Sequence contig(branch);

	size_t kmerCount = branch.calculateBranchMultiplicity();
	if (writer != NULL)
		writer->WriteSequence(contig, id, kmerCount);

	// Remove low-coverage contigs.
	float coverage = (float)kmerCount / branch.size();
	if (opt::coverage > 0 && coverage < opt::coverage) {
		for (BranchRecord::iterator it = branch.begin();
				it != branch.end(); ++it)
			put(vertex_removed, g, it->first, true);
		return branch.size();
	}
	return 0;
}

/** Assemble contigs.
 * @return the number of contigs assembled
 */
template<typename Graph>
size_t
assemble(Graph& g, FastaWriter* fileWriter)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::vertex_iterator Vit;
	typedef typename vertex_bundle_type<Graph>::type VP;

	Timer timer("Assemble");

	size_t kmerCount = 0;
	unsigned contigID = 0;
	size_t assembledKmer = 0;
	size_t lowCoverageKmer = 0;
	size_t lowCoverageContigs = 0;

	std::pair<Vit, Vit> urange = vertices(g);
	for (Vit uit = urange.first; uit != urange.second; ++uit) {
		V u = *uit;
		const VP& up = get(vertex_bundle, g, u);
		if (up.deleted())
			continue;

		if (++kmerCount % 1000000 == 0)
			logger(1) << "Assembling: " << kmerCount << '\n';

		extDirection dir;
		SeqContiguity status = checkSeqContiguity(up, dir, true);
		if (status == SC_CONTIGUOUS)
			continue;
		else if(status == SC_ISLAND)
		{
			BranchRecord currBranch(SENSE);
			currBranch.push_back(std::make_pair(u, up));
			currBranch.terminate(BS_NOEXT);
			size_t removed = assembleContig(g,
					fileWriter, currBranch, contigID++);
			assembledKmer += currBranch.size();
			if (removed > 0) {
				lowCoverageContigs++;
				lowCoverageKmer += removed;
			}
			continue;
		}
		assert(status == SC_ENDPOINT);

		BranchRecord currBranch(dir);
		currBranch.push_back(std::make_pair(u, up));
		Kmer v = u;
		extendBranch(currBranch, v, up.getExtension(dir));
		assert(currBranch.isActive());
		while (currBranch.isActive()) {
			const VP& vp = g[orientVertex(g, v)];
			processLinearExtensionForBranch(currBranch, v,
					vp.extension(), vp.getMultiplicity(),
					UINT_MAX);
		}

		if (currBranch.isCanonical()) {
			size_t removed = assembleContig(g,
					fileWriter, currBranch, contigID++);
			assembledKmer += currBranch.size();
			if (removed > 0) {
				lowCoverageContigs++;
				lowCoverageKmer += removed;
			}
		}

		pumpNetwork(g);
	}

	if (opt::coverage > 0) {
		std::cout << "Found " << assembledKmer << " k-mer in "
			<< contigID
			<< " contigs before removing low-coverage contigs.\n"
			"Removed " << lowCoverageKmer << " k-mer in "
				<< lowCoverageContigs << " low-coverage contigs.\n";
	} else {
		assert(assembledKmer <= kmerCount);
		size_t circularKmer = kmerCount - assembledKmer;
		if (circularKmer > 0)
			std::cout << "Left " << circularKmer
				<< " unassembled k-mer in circular contigs.\n";
		std::cout << "Assembled " << assembledKmer << " k-mer in "
			<< contigID << " contigs.\n";
	}
	return contigID;
}

} // namespace dbg

#endif
