#ifndef ASSEMBLY_SPLITALGORITHM_H
#define ASSEMBLY_SPLITALGORITHM_H 1

namespace AssemblyAlgorithms {

/** Mark the specified vertex and its neighbours.
 * @return the number of marked edges
 */
template <typename Graph>
size_t markNeighbours(Graph* g,
		const typename Graph::value_type& u, extDirection sense)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename std::vector<V> Vector;

	Vector adj;
	generateSequencesFromExtension(u.first, sense,
			u.second.getExtension(sense), adj);
	for (typename Vector::iterator v = adj.begin(); v != adj.end(); ++v)
		g->mark(*v, !sense);
	return adj.size();
}

/** Mark ambiguous branches and branches from palindromes for removal.
 * @return the number of branches marked
 */
template <typename Graph>
size_t markAmbiguous(Graph* g)
{
	typedef typename Graph::iterator iterator;

	Timer timer(__func__);
	size_t progress = 0;
	size_t countv = 0, counte = 0;
	for (iterator it = g->begin(); it != g->end(); ++it) {
		if (it->second.deleted())
			continue;

		if (++progress % 1000000 == 0)
			logger(1) << "Splitting: " << progress << '\n';

		if (!opt::ss && it->first.isPalindrome()) {
			countv += 2;
			g->mark(it->first);
			counte += markNeighbours(g, *it, SENSE);
		} else {
			for (extDirection sense = SENSE;
					sense <= ANTISENSE; ++sense) {
				if (it->second.getExtension(sense).isAmbiguous()
						|| (!opt::ss && it->first.isPalindrome(sense))) {
					countv++;
					g->mark(it->first, sense);
					counte += markNeighbours(g, *it, sense);
				}
			}
		}

		g->pumpNetwork();
	}
	tempCounter[5] = countv;
	logger(0) << "Marked " << counte << " edges of " << countv
		<< " ambiguous vertices." << std::endl;

	return countv;
}

/** Remove the edges of marked and deleted vertices.
 * @return the number of branches removed
 */
template <typename Graph>
size_t splitAmbiguous(Graph* pSC)
{
	typedef typename Graph::iterator iterator;

	Timer timer(__func__);
	size_t count = 0;
	for (iterator it = pSC->begin();
			it != pSC->end(); ++it) {
		if (!it->second.deleted())
			continue;
		for (extDirection sense = SENSE;
				sense <= ANTISENSE; ++sense) {
			if (it->second.marked(sense)) {
				removeExtensionsToSequence(pSC, *it, sense);
				count++;
			}
		}
		pSC->pumpNetwork();
	}
	tempCounter[7] += count;
	logger(0) << "Split " << count << " ambigiuous branches.\n";
	return count;
}

} // namespace AssemblyAlgorithms

#endif
