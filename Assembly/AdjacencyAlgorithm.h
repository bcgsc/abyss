#ifndef ASSEMBLY_ADJACENCYALGORITHM_H
#define ASSEMBLY_ADJACENCYALGORITHM_H 1

namespace AssemblyAlgorithms {

/** Generate the adjacency information for each sequence in the
 * collection. */
template <typename Graph>
size_t generateAdjacency(Graph* seqCollection)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename Graph::Symbol Symbol;
	typedef typename Graph::SymbolSet SymbolSet;

	Timer timer("GenerateAdjacency");

	size_t count = 0;
	size_t numBasesSet = 0;
	for (typename Graph::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		if (++count % 1000000 == 0)
			logger(1) << "Finding adjacent k-mer: " << count << '\n';

		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			V testSeq(iter->first);
			Symbol adjBase = testSeq.shift(dir);
			for (unsigned i = 0; i < SymbolSet::NUM; ++i) {
				testSeq.setLastBase(dir, Symbol(i));
				if (seqCollection->setBaseExtension(
							testSeq, !dir, adjBase))
					numBasesSet++;
			}
		}
		seqCollection->pumpNetwork();
	}

	if (numBasesSet > 0) {
		logger(0) << "Added " << numBasesSet << " edges.\n";
		if (!opt::db.empty())
			addToDb("EdgesGenerated", numBasesSet);
	}
	return numBasesSet;
}

} // namespace AssemblyAlgorithms

#endif
