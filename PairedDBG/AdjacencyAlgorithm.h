#ifndef PAIREDDBG_ADJACENCYALGORITHM_H
#define PAIREDDBG_ADJACENCYALGORITHM_H 1

namespace AssemblyAlgorithms {

/** Generate the adjacency information for each sequence in the
 * collection. */
static inline
size_t generateAdjacency(ISequenceCollection* seqCollection)
{
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;

	Timer timer("GenerateAdjacency");

	size_t count = 0;
	size_t numBasesSet = 0;
	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		if (++count % 1000000 == 0)
			logger(1) << "Finding adjacent k-mer: " << count << '\n';

		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			V testSeq(iter->first);
			Dinuc adjBase = testSeq.shift(dir);
			for (unsigned i = 0; i < DinucSet::NUM_EDGES; i++) {
				testSeq.setLastBase(dir, i);
				if (seqCollection->setBaseExtension(
							testSeq, !dir, adjBase))
					numBasesSet++;
			}
		}
		seqCollection->pumpNetwork();
	}

	if (numBasesSet > 0) {
		logger(0) << "Added " << numBasesSet << " edges.\n";
#if _SQL
		addToDb("EdgesGenerated", numBasesSet);
#endif
	}
	return numBasesSet;
}

} // namespace AssemblyAlgorithms

#endif
