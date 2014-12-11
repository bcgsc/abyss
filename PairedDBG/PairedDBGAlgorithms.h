#ifndef PairedDBG_PairedDBGAlgorithms_h
#define PairedDBG_PairedDBGAlgorithms_h 1

#include "Graph/GraphAlgorithms.h" // for removeEdgeIf
#include <iostream>

/** Return true if a paired dBG edge is inconsistent.
 * This predicate is only valid when the size of the gap is exactly zero.
 */
template <typename Graph>
struct InconsistentPairedDBGEdge
{
	InconsistentPairedDBGEdge(Graph& g) : m_g(g) { }
	bool operator()(typename graph_traits<Graph>::edge_descriptor e) const
	{
		assert(opt::kmerSize == 2 * opt::singleKmerSize);
		// u aaaaabbbbb
		// v  aaaaabbbbb
		return source(e, m_g).front().b()
			!= target(e, m_g).back().a();
	}
	const Graph& m_g;
};

/** Remove inconsistent edges from a paired dBG.
 * Currently this algorithm is only used for the special case when
 * assembling a paired dBG whose gap is exactly zero.
 */
template <typename Graph>
void
removePairedDBGInconsistentEdges(Graph& g)
{
	assert(opt::kmerSize >= 2 * opt::singleKmerSize);
	if (opt::kmerSize == 2 * opt::singleKmerSize) {
		// The gap is exactly zero.
		std::cerr << "Removed "
			<< removeEdgeIf(InconsistentPairedDBGEdge<Graph>(g), g)
			<< " inconsistent edges.\n";
	}
}

#endif
