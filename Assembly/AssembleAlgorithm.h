#ifndef ASSEMBLY_ASSEMBLEALGORITHM_H
#define ASSEMBLY_ASSEMBLEALGORITHM_H 1

#include "DataLayer/FastaWriter.h"
#include <iostream>
#include <climits>

namespace AssemblyAlgorithms {

/** Assemble a contig.
 * @return the number of k-mer below the coverage threshold
 */
template <typename Graph>
size_t assembleContig(
		Graph* seqCollection, FastaWriter* writer,
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
			seqCollection->remove(it->first);
		return branch.size();
	}
	return 0;
}

/** Assemble contigs.
 * @return the number of contigs assembled
 */
static inline
size_t assemble(SequenceCollectionHash* seqCollection,
		FastaWriter* fileWriter = NULL)
{
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef Graph::SymbolSetPair SymbolSetPair;

	Timer timer("Assemble");

	size_t kmerCount = 0;
	unsigned contigID = 0;
	size_t assembledKmer = 0;
	size_t lowCoverageKmer = 0;
	size_t lowCoverageContigs = 0;

	for (Graph::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;
		kmerCount++;

		extDirection dir;
		SeqContiguity status = checkSeqContiguity(*iter, dir, true);
		if (status == SC_CONTIGUOUS)
			continue;
		else if(status == SC_ISLAND)
		{
			BranchRecord currBranch(SENSE);
			currBranch.push_back(*iter);
			currBranch.terminate(BS_NOEXT);
			size_t removed = assembleContig(seqCollection,
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
		currBranch.push_back(*iter);
		V currSeq = iter->first;
		extendBranch(currBranch, currSeq,
				iter->second.getExtension(dir));
		assert(currBranch.isActive());
		while(currBranch.isActive())
		{
			SymbolSetPair extRec;
			int multiplicity = -1;
			bool success = seqCollection->getSeqData(
					currSeq, extRec, multiplicity);
			assert(success);
			(void)success;
			processLinearExtensionForBranch(currBranch,
					currSeq, extRec, multiplicity, UINT_MAX);
		}

		if ((opt::ss && currBranch.getDirection() == SENSE)
				|| (!opt::ss && currBranch.isCanonical())) {
			size_t removed = assembleContig(seqCollection,
					fileWriter, currBranch, contigID++);
			assembledKmer += currBranch.size();
			if (removed > 0) {
				lowCoverageContigs++;
				lowCoverageKmer += removed;
			}
		}

		seqCollection->pumpNetwork();
	}

	if (opt::coverage > 0) {
		std::cout << "Found " << assembledKmer << " k-mer in " << contigID
			<< " contigs before removing low-coverage contigs.\n"
			"Removed " << lowCoverageKmer << " k-mer in "
				<< lowCoverageContigs << " low-coverage contigs.\n";
		tempCounter[3] += lowCoverageContigs;
		tempCounter[4] += lowCoverageKmer;
	} else {
		assert(assembledKmer <= kmerCount);
		size_t circularKmer = kmerCount - assembledKmer;
		if (circularKmer > 0)
			std::cout << "Left " << circularKmer
				<< " unassembled k-mer in circular contigs.\n";
		std::cout << "Assembled " << assembledKmer << " k-mer in "
			<< contigID << " contigs.\n";
		if (!opt::db.empty()) {
			addToDb("finalAmbgVertices", tempCounter[5]);
			//addToDb("finalAmbgEdges", tempCounter[6]);
			tempCounter.assign(16,0);
			addToDb("assembledKmerNum", assembledKmer);
			addToDb("assembledCntg", contigID);
		}
	}
	return contigID;
}

} // namespace AssemblyAlgorithms

#endif
