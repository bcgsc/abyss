#ifndef ASSEMBLY_BUBBLEALGORITHM_H
#define ASSEMBLY_BUBBLEALGORITHM_H 1

#include "Common/IOUtil.h"
#include "Common/Options.h" // for opt::rank
#include <fstream>
#include <sstream>

namespace AssemblyAlgorithms {

static inline bool
processBranchGroupExtension(BranchGroup& group,
		size_t branchIndex,
		const graph_traits<SequenceCollectionHash>::vertex_descriptor& seq,
		SequenceCollectionHash::SymbolSetPair ext,
		int multiplicity,
		unsigned maxLength);

template <typename Graph>
void collapseJoinedBranches(Graph* collection,
		BranchGroup& group);

static inline
void writeBubble(std::ostream& out, const BranchGroup& group, unsigned id);

/** Open the bubble file. */
static inline
void openBubbleFile(std::ofstream& out)
{
	if (opt::snpPath.empty())
		return;
	std::string path;
	if (opt::rank < 0) {
		path = opt::snpPath;
	} else {
		std::ostringstream s;
		s << "snp-" << opt::rank << ".fa";
		path = s.str();
	}
	out.open(path.c_str());
	assert_good(out, path);
}

/** Pop bubbles. */
static inline
size_t popBubbles(SequenceCollectionHash* seqCollection, std::ostream& out)
{
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef Graph::SymbolSetPair SymbolSetPair;

	Timer timer("PopBubbles");
	size_t numPopped = 0;

	// Set the cutoffs
	const unsigned maxNumBranches = 3;
	const unsigned maxLength = opt::bubbleLen - V::length() + 1;

	for (Graph::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		SymbolSetPair extRec = iter->second.extension();
		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			if (extRec.dir[dir].isAmbiguous()) {
				// Found a potential bubble, examine each branch
				bool stop = false;

				// Create the branch group
				BranchGroup branchGroup(dir, maxNumBranches,
						iter->first);
				initiateBranchGroup(branchGroup, iter->first,
						extRec.dir[dir]);

				// Iterate over the branches
				while(!stop)
				{
					size_t numBranches = branchGroup.size();
					for (unsigned j = 0; j < numBranches; ++j) {
						// Get the extensions of this branch
						SymbolSetPair extRec;
						int multiplicity = -1;

						const V& lastKmer
							= branchGroup[j].back().first;
						bool success = seqCollection->getSeqData(
								lastKmer, extRec, multiplicity);
						assert(success);
						(void)success;
						processBranchGroupExtension(branchGroup, j,
								lastKmer, extRec, multiplicity,
								maxLength);
					}

					// At this point all branches should have the same
					// length or one will be a noext.
					branchGroup.updateStatus(maxLength);
					BranchGroupStatus status
						= branchGroup.getStatus();
					if (status == BGS_TOOLONG
							|| status == BGS_TOOMANYBRANCHES
							|| status == BGS_NOEXT) {
						stop = true;
					}
					else if(status == BGS_JOINED)
					{
						static unsigned snpID;
						writeBubble(out, branchGroup, ++snpID);
						assert(branchGroup.isAmbiguous(
									*seqCollection));
						collapseJoinedBranches(seqCollection,
								branchGroup);
						assert(!branchGroup.isAmbiguous(
									*seqCollection));
						numPopped++;
						stop = true;
					} else
						assert(status == BGS_ACTIVE);
				}
			}
		}
		seqCollection->pumpNetwork();
	}

	if (numPopped > 0)
		std::cout << "Removed " << numPopped << " bubbles.\n";
	if (!opt::db.empty()) {
		addToDb("totalErodedTips", tempCounter[0]);
		addToDb("totalPrunedTips", tempCounter[1]);
		addToDb("totalLowCovCntg", tempCounter[3]);
		addToDb("totalLowCovKmer", tempCounter[4]);
		addToDb("totalSplitAmbg", tempCounter[7]);
		addToDb("poppedBubbles", numPopped);
		tempCounter.assign(16,0);
	}
	return numPopped;
}

// Populate a branch group with the inital branches from a sequence
static inline void
initiateBranchGroup(BranchGroup& group,
		const graph_traits<SequenceCollectionHash>::vertex_descriptor& seq,
		const SequenceCollectionHash::SymbolSet& extension)
{
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;

	std::vector<V> extSeqs;
	generateSequencesFromExtension(seq, group.getDirection(),
			extension, extSeqs);
	assert(extSeqs.size() > 1);
	for (std::vector<V>::iterator seqIter = extSeqs.begin();
			seqIter != extSeqs.end(); ++seqIter)
		group.addBranch(BranchRecord(group.getDirection()), *seqIter);
}

/** Process an a branch group extension. */
static inline bool
processBranchGroupExtension(BranchGroup& group,
		size_t branchIndex,
		const graph_traits<SequenceCollectionHash>::vertex_descriptor& seq,
		SequenceCollectionHash::SymbolSetPair ext,
		int multiplicity,
		unsigned maxLength)
{
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef vertex_bundle_type<Graph>::type VP;

	BranchRecord& branch = group[branchIndex];
	branch.setData(std::make_pair(seq, VP(multiplicity, ext)));

	extDirection dir = group.getDirection();
	if (ext.dir[!dir].isAmbiguous()) {
		// Check that this fork is due to branches of our bubble
		// merging back together. If not, stop this bubble.
		if (branch.size() < 2) {
			group.setNoExtension();
			return false;
		}

		std::vector<V> extKmer;
		generateSequencesFromExtension(seq, !dir,
				ext.dir[!dir], extKmer);
		assert(extKmer.size() > 1);
		for (std::vector<V>::iterator it = extKmer.begin();
				it != extKmer.end(); ++it) {
			assert(branch.size() > 1);
			if (!group.exists(branch.size() - 2, *it)) {
				group.setNoExtension();
				return false;
			}
		}
		// Ignore the ambiguity.
		ext.dir[!dir].clear();
	}

	if (ext.dir[dir].isAmbiguous()) {
		// Create a new branch to follow the fork.
		std::vector<V> extKmer;
		generateSequencesFromExtension(seq, dir,
				ext.dir[dir], extKmer);
		assert(extKmer.size() > 1);
		BranchRecord original = branch;
		std::vector<V>::iterator it = extKmer.begin();
		branch.push_back(std::make_pair(*it++, VP()));
		for (; it != extKmer.end(); ++it)
			group.addBranch(original, *it);
		return group.isExtendable();
	}

	V nextKmer = seq;
	if (processLinearExtensionForBranch(branch,
			nextKmer, ext, multiplicity,
			maxLength, false))
		branch.push_back(std::make_pair(nextKmer, VP()));
	else
		group.setNoExtension();
	return group.isExtendable();
}

/** Write a bubble to the specified file. */
static inline
void writeBubble(std::ostream& out, const BranchGroup& group, unsigned id)
{
	if (opt::snpPath.empty())
		return;

	char allele = 'A';
	for (BranchGroup::const_iterator it = group.begin();
			it != group.end(); ++it) {
		const BranchRecord& currBranch = *it;
		Sequence contig(currBranch);
		out << '>' << id << allele++ << ' '
			<< contig.length() << ' '
			<< currBranch.calculateBranchMultiplicity() << '\n'
			<< contig.c_str() << '\n';
	}
	assert_good(out, opt::snpPath);
}

/** Collapse a bubble to a single path. */
template <typename Graph>
void collapseJoinedBranches(Graph* collection,
		BranchGroup& group)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename vertex_bundle_type<Graph>::type VP;
	typedef typename std::map<V, VP> Map;

	const BranchRecord& best = group[0];
	logger(5) << "Popping " << best.size() << ' '
		<< best.front().first << '\n';

	// Add the k-mer from the dead branches.
	Map doomed;
	for (BranchGroup::const_iterator branchIt = group.begin() + 1;
			branchIt != group.end(); ++branchIt) {
		const BranchRecord& branch = *branchIt;
		for (BranchRecord::const_iterator it = branch.begin();
				it != branch.end(); ++it)
			doomed.insert(*it);
	}

	// Remove the k-mer that are in the good branch.
	for (BranchRecord::const_iterator it = best.begin();
			it != best.end(); ++it)
		doomed.erase(it->first);

	// Remove the dead k-mer from the assembly.
	for (typename Map::const_iterator it = doomed.begin();
			it != doomed.end(); ++it)
		removeSequenceAndExtensions(collection, *it);
}

} // namespace AssemblyAlgorithms

#endif
