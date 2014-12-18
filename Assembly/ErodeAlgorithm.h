#ifndef ASSEMBLY_ERODE_ALGORITHM
#define ASSEMBLY_ERODE_ALGORITHM 1

namespace AssemblyAlgorithms {

/** The number of k-mer that have been eroded. */
extern size_t g_numEroded;

template <typename Graph>
void removeExtensionsToSequence(Graph* seqCollection,
		const typename Graph::value_type& seq, extDirection dir);

/**
 * Remove a k-mer and update the extension records of the k-mer that
 * extend to it.
 */
template <typename Graph>
void removeSequenceAndExtensions(Graph* seqCollection,
		const typename Graph::value_type& seq)
{
	// This removes the reverse complement as well
	seqCollection->remove(seq.first);
	removeExtensionsToSequence(seqCollection, seq, SENSE);
	removeExtensionsToSequence(seqCollection, seq, ANTISENSE);
}

/** Remove all the extensions to this sequence. */
template <typename Graph>
void removeExtensionsToSequence(Graph* seqCollection,
		const typename Graph::value_type& seq, extDirection dir)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename Graph::Symbol Symbol;
	typedef typename Graph::SymbolSet SymbolSet;

	SymbolSet extension(seq.second.getExtension(dir));
	V testSeq(seq.first);
	Symbol extBase = testSeq.shift(dir);
	for (unsigned i = 0; i < extension.NUM; ++i) {
		Symbol x(i);
		if (extension.checkBase(x)) {
			testSeq.setLastBase(dir, x);
			seqCollection->removeExtension(testSeq, !dir, extBase);
		}
	}
}

/** Return the number of k-mer that have been eroded. */
static inline
size_t getNumEroded()
{
	size_t numEroded = g_numEroded;
	g_numEroded = 0;
	tempCounter[0] += numEroded;
	logger(0) << "Eroded " << numEroded << " tips.\n";
	return numEroded;
}

/** Consider the specified k-mer for erosion.
 * @return the number of k-mer eroded, zero or one
 */
template <typename Graph>
size_t erode(Graph* c, const typename Graph::value_type& seq)
{
	typedef typename vertex_bundle_type<Graph>::type VP;

	if (seq.second.deleted())
		return 0;
	extDirection dir;
	SeqContiguity contiguity = checkSeqContiguity(seq, dir);
	if (contiguity == SC_CONTIGUOUS)
		return 0;

	const VP& data = seq.second;
	if (data.getMultiplicity() < opt::erode
			|| data.getMultiplicity(SENSE) < opt::erodeStrand
			|| data.getMultiplicity(ANTISENSE) < opt::erodeStrand) {
		removeSequenceAndExtensions(c, seq);
		g_numEroded++;
		return 1;
	} else
		return 0;
}

/** The given sequence has changed. */
static inline
void erosionObserver(SequenceCollectionHash* c,
		const SequenceCollectionHash::value_type& seq)
{
	erode(c, seq);
}

//
// Erode data off the ends of the graph, one by one
//
template <typename Graph>
size_t erodeEnds(Graph* seqCollection)
{
	typedef typename Graph::iterator iterator;

	Timer erodeEndsTimer("Erode");
	assert(g_numEroded == 0);
	seqCollection->attach(erosionObserver);

	for (iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		erode(seqCollection, *iter);
		seqCollection->pumpNetwork();
	}

	seqCollection->detach(erosionObserver);
	return getNumEroded();
}

} // namespace AssemblyAlgorithms

#endif
