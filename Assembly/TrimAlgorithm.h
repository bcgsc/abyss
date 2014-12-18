#ifndef ASSEMBLY_TRIMALGORITHM_H
#define ASSEMBLY_TRIMALGORITHM_H 1

namespace AssemblyAlgorithms {

template <typename Graph>
bool processTerminatedBranchTrim(Graph* seqCollection, BranchRecord& branch);

static inline
size_t trimSequences(SequenceCollectionHash* seqCollection,
		unsigned maxBranchCull);

/** Trimming driver function */
static inline
void performTrim(SequenceCollectionHash* seqCollection)
{
	if (opt::trimLen == 0)
		return;
	unsigned rounds = 0;
	size_t total = 0;
	for (unsigned trim = 1; trim < opt::trimLen; trim *= 2) {
		rounds++;
		total += trimSequences(seqCollection, trim);
	}
	size_t count;
	while ((count = trimSequences(seqCollection, opt::trimLen)) > 0) {
		rounds++;
		total += count;
	}
	std::cout << "Pruned " << total << " tips in "
		<< rounds << " rounds.\n";
	tempCounter[1] += total;
	tempCounter[2] = rounds;
}

/** Prune tips shorter than maxBranchCull. */
static inline
size_t trimSequences(SequenceCollectionHash* seqCollection,
		unsigned maxBranchCull)
{
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef Graph::SymbolSetPair SymbolSetPair;

	Timer timer("TrimSequences");
	std::cout << "Pruning tips shorter than "
		<< maxBranchCull << " bp...\n";
	size_t numBranchesRemoved = 0;

	for (Graph::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		extDirection dir;
		// dir will be set to the trimming direction if the sequence
		// can be trimmed.
		SeqContiguity status = checkSeqContiguity(*iter, dir);

		if (status == SC_CONTIGUOUS)
			continue;
		else if(status == SC_ISLAND)
		{
			// remove this sequence, it has no extensions
			seqCollection->mark(iter->first);
			numBranchesRemoved++;
			continue;
		}

		BranchRecord currBranch(dir);
		V currSeq = iter->first;
		while(currBranch.isActive())
		{
			SymbolSetPair extRec;
			int multiplicity = -1;
			bool success = seqCollection->getSeqData(
					currSeq, extRec, multiplicity);
			assert(success);
			(void)success;
			processLinearExtensionForBranch(currBranch,
					currSeq, extRec, multiplicity, maxBranchCull);
		}

		// The branch has ended check it for removal, returns true if
		// it was removed.
		if(processTerminatedBranchTrim(seqCollection, currBranch))
		{
			numBranchesRemoved++;
		}
		seqCollection->pumpNetwork();
	}

	size_t numSweeped = removeMarked(seqCollection);

	if (numBranchesRemoved > 0)
		logger(0) << "Pruned " << numSweeped << " k-mer in "
			<< numBranchesRemoved << " tips.\n";
	return numBranchesRemoved;
}

/** Extend this branch. */
static inline
bool extendBranch(BranchRecord& branch,
		graph_traits<SequenceCollectionHash>::vertex_descriptor& kmer,
		SequenceCollectionHash::SymbolSet ext)
{
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;
	
	if (!ext.hasExtension()) {
		branch.terminate(BS_NOEXT);
		return false;
	} else if (ext.isAmbiguous()) {
		branch.terminate(BS_AMBI_SAME);
		return false;
	} else {
		std::vector<V> adj;
		generateSequencesFromExtension(kmer, branch.getDirection(),
				ext, adj);
		assert(adj.size() == 1);
		kmer = adj.front();
		return true;
	}
}

/**
 * Process the extension for this branch for the trimming algorithm
 * CurrSeq is the current sequence being inspected (the next member to
 * be added to the branch). The extension record is the extensions of
 * that sequence and multiplicity is the number of times that kmer
 * appears in the data set. After processing currSeq is unchanged if
 * the branch is no longer active or else it is the generated
 * extension. If the parameter addKmer is true, add the k-mer to the
 * branch.
 */
static inline bool
processLinearExtensionForBranch(BranchRecord& branch,
		graph_traits<SequenceCollectionHash>::vertex_descriptor& currSeq,
		SequenceCollectionHash::SymbolSetPair extensions,
		int multiplicity,
		unsigned maxLength, bool addKmer)
{
	typedef SequenceCollectionHash Graph;
	typedef vertex_bundle_type<Graph>::type VP;

	/** Stop contig assembly at palindromes. */
	const bool stopAtPalindromes = !opt::ss && maxLength == UINT_MAX;

	extDirection dir = branch.getDirection();
	if (branch.isTooLong(maxLength)) {
		// Too long.
		branch.terminate(BS_TOO_LONG);
		return false;
	} else if (extensions.dir[!dir].isAmbiguous()) {
		// Ambiguous.
		branch.terminate(BS_AMBI_OPP);
		return false;
	} else if (stopAtPalindromes && currSeq.isPalindrome()) {
		// Palindrome.
		branch.terminate(BS_AMBI_SAME);
		return false;
	}

	if (addKmer)
		branch.push_back(std::make_pair(currSeq,
					VP(multiplicity, extensions)));

	if (branch.isTooLong(maxLength)) {
		// Too long.
		branch.terminate(BS_TOO_LONG);
		return false;
	} else if (stopAtPalindromes && currSeq.isPalindrome(dir)) {
		// Palindrome.
		branch.terminate(BS_AMBI_SAME);
		return false;
	}

	return extendBranch(branch, currSeq, extensions.dir[dir]);
}

/** Trim the specified branch if it meets trimming criteria.
 * @return true if the specified branch was trimmed
 */
template <typename Graph>
bool processTerminatedBranchTrim(Graph* seqCollection, BranchRecord& branch)
{
	assert(!branch.isActive());
	assert(!branch.empty());
	if (branch.getState() == BS_NOEXT
			|| branch.getState() == BS_AMBI_OPP) {
		logger(5) << "Pruning " << branch.size() << ' '
			<< branch.front().first << '\n';
		for (BranchRecord::iterator it = branch.begin();
				it != branch.end(); ++it)
			seqCollection->mark(it->first);
		return true;
	} else
		return false;
}

} // namespace AssemblyAlgorithms

#endif
