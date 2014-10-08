#ifndef PAIREDDBG_ASSEMBLEALGORITHM_H
#define PAIREDDBG_ASSEMBLEALGORITHM_H 1

#include <iostream>

namespace AssemblyAlgorithms {

/** Assemble a contig.
 * @return the number of k-mer below the coverage threshold
 */
static inline
size_t assembleContig(
		ISequenceCollection* seqCollection, FastaWriter* writer,
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
		FastaWriter* fileWriter)
{
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;

	Timer timer("Assemble");

	size_t kmerCount = 0;
	unsigned contigID = 0;
	size_t assembledKmer = 0;
	size_t lowCoverageKmer = 0;
	size_t lowCoverageContigs = 0;

	for (ISequenceCollection::iterator iter = seqCollection->begin();
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
			DinucSetPair extRec;
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
#if _SQL
		tempCounter[3] += lowCoverageContigs;
		tempCounter[4] += lowCoverageKmer;
#endif
	} else {
		assert(assembledKmer <= kmerCount);
		size_t circularKmer = kmerCount - assembledKmer;
		if (circularKmer > 0)
			std::cout << "Left " << circularKmer
				<< " unassembled k-mer in circular contigs.\n";
		std::cout << "Assembled " << assembledKmer << " k-mer in "
			<< contigID << " contigs.\n";
#if _SQL
		addToDb("finalAmbgVertices", tempCounter[5]);
		//addToDb("finalAmbgEdges", tempCounter[6]);
		tempCounter.assign(16,0);
		addToDb("assembledKmerNum", assembledKmer);
		addToDb("assembledCntg", contigID);
#endif
	}
	return contigID;
}

} // namespace AssemblyAlgorithms

#endif
