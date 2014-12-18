#ifndef ASSEMBLY_ASSEMBLYALGORITHMS_H
#define ASSEMBLY_ASSEMBLYALGORITHMS_H 1

#include "Assembly/BranchGroup.h"
#include "Assembly/Options.h"
#include "Common/Log.h"
#include "Common/Timer.h"
#include <vector>
#include <string>
#include "Common/InsOrderedMap.h"

class Histogram;

/** A summary of the in- and out-degree of a vertex. */
enum SeqContiguity
{
	SC_ISLAND, // sequence is completely isolated
	SC_ENDPOINT, // one end of the sequence is open
	SC_CONTIGUOUS // the sequence is closed on both ends
};

/** De Bruijn graph assembly algorithms. */
namespace AssemblyAlgorithms {

extern std::vector<size_t> tempCounter;
extern InsOrderedMap<std::string,int> tempStatMap;
extern void addToDb(const std::string&, const int&);

static inline
bool extendBranch(BranchRecord& branch,
		graph_traits<SequenceCollectionHash>::vertex_descriptor& kmer,
		SequenceCollectionHash::SymbolSet ext);

static inline bool
processLinearExtensionForBranch(BranchRecord& branch,
		graph_traits<SequenceCollectionHash>::vertex_descriptor& currSeq,
		SequenceCollectionHash::SymbolSetPair extensions,
		int multiplicity,
		unsigned maxLength, bool addKmer = true);

static inline void
initiateBranchGroup(BranchGroup& group,
		const graph_traits<SequenceCollectionHash>::vertex_descriptor& seq,
		const SequenceCollectionHash::SymbolSet& extension);

template <typename Graph>
void removeSequenceAndExtensions(Graph* seqCollection,
		const typename Graph::value_type& seq);

/** Return the kmer which are adjacent to this kmer. */
template <typename V, typename SymbolSet>
void generateSequencesFromExtension(
		const V& currSeq,
		extDirection dir,
		SymbolSet extension,
		std::vector<V>& outseqs)
{
	typedef typename SymbolSet::Symbol Symbol;

	std::vector<V> extensions;
	V extSeq(currSeq);
	extSeq.shift(dir);

	// Check for the existance of the 4 possible extensions
	for (unsigned i = 0; i < SymbolSet::NUM; ++i) {
		// Does this sequence have an extension?
		Symbol x(i);
		if (extension.checkBase(x)) {
			extSeq.setLastBase(dir, x);
			outseqs.push_back(extSeq);
		}
	}
}

/** Return the adjacency of this sequence.
 * @param considerMarks when true, treat a marked vertex as having
 * no edges
 */
static inline
SeqContiguity checkSeqContiguity(
		const SequenceCollectionHash::value_type& seq,
		extDirection& outDir, bool considerMarks = false)
{
	assert(!seq.second.deleted());
	bool child = seq.second.hasExtension(SENSE)
		&& !(considerMarks && seq.second.marked(SENSE));
	bool parent = seq.second.hasExtension(ANTISENSE)
		&& !(considerMarks && seq.second.marked(ANTISENSE));
	if(!child && !parent)
	{
		//this sequence is completely isolated
		return SC_ISLAND;
	}
	else if(!child)
	{
		outDir = ANTISENSE;
		return SC_ENDPOINT;
	}
	else if(!parent)
	{
		outDir = SENSE;
		return SC_ENDPOINT;
	}
	else
	{
		// sequence is contiguous
		return SC_CONTIGUOUS;
	}
}

/** Remove all marked k-mer.
 * @return the number of removed k-mer
 */
template <typename Graph>
size_t removeMarked(Graph* pSC)
{
	typedef typename Graph::iterator iterator;

	Timer timer(__func__);
	size_t count = 0;
	for (iterator it = pSC->begin(); it != pSC->end(); ++it) {
		if (it->second.deleted())
			continue;
		if (it->second.marked()) {
			removeSequenceAndExtensions(pSC, *it);
			count++;
		}
		pSC->pumpNetwork();
	}
	if (count > 0)
		logger(1) << "Removed " << count << " marked k-mer.\n";
	return count;
}

} // namespace AssemblyAlgorithms

#include "AdjacencyAlgorithm.h"
#include "AssembleAlgorithm.h"
#include "BubbleAlgorithm.h"
#include "CoverageAlgorithm.h"
#include "ErodeAlgorithm.h"
#include "LoadAlgorithm.h"
#include "SplitAlgorithm.h"
#include "TrimAlgorithm.h"

#endif
