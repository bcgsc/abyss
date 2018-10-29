#ifndef _ROLLING_BLOOM_DBG_VISITOR_H_
#define _ROLLING_BLOOM_DBG_VISITOR_H_ 1

#include "Common/UnorderedMap.h"
#include "Graph/BreadthFirstSearch.h"

#include <iostream>

/**
 * Visitor class that outputs visited nodes/edges in GraphViz format during
 * a bounded breadth first traversal. An instance of this class may be passed
 * as an argument to the `breadthFirstSearch` function.
 */
template <typename GraphT>
class RollingBloomDBGVisitor : public DefaultBFSVisitor<GraphT>
{
public:

	typedef typename boost::graph_traits<GraphT>::vertex_descriptor VertexT;
	typedef typename boost::graph_traits<GraphT>::edge_descriptor EdgeT;

	typedef unordered_map<VertexT, unsigned> DepthMap;
	typedef typename DepthMap::const_iterator DepthMapConstIt;
	typedef typename DepthMap::iterator DepthMapIt;

	typedef std::vector<std::pair<std::string, unordered_set<VertexT> > >
		KmerProperties;
	typedef typename KmerProperties::const_iterator KmerPropertiesIt;

	/** Constructor */
	RollingBloomDBGVisitor(const VertexT& root, unsigned maxDepth,
		const KmerProperties& kmerProperties, std::ostream& out)
		: m_root(root), m_maxDepth(maxDepth),
		m_kmerProperties(kmerProperties), m_out(out)
	{
		m_depthMap[root] = 0;

		/* start directed graph (GraphViz) */
		m_out << "digraph " << root.kmer().c_str() << " {\n";
	}

	/** Called when graph search completes */
	void post_processing()
	{
		/* end directed graph (GraphViz) */
		m_out << "}\n";
	}

	/** Invoked on edges that lead to an undiscovered vertex */
	BFSVisitorResult tree_edge(const EdgeT& e, const GraphT& g)
	{
		VertexT parent = source(e, g);
		VertexT child = target(e, g);

		if (!isForwardEdge(e, g))
			std::swap(parent, child);

		DepthMapConstIt parentIt = m_depthMap.find(parent);
		assert(parentIt != m_depthMap.end());
		unsigned parentDepth = parentIt->second;

		if (parentDepth >= m_maxDepth)
			return BFS_SKIP_ELEMENT;

		DepthMapIt childIt;
		bool inserted;
		boost::tie(childIt, inserted) = m_depthMap.insert(
			std::make_pair(child, parentDepth + 1));
		assert(inserted);

		/*
		 * since we use this visitor class with an undirected BFS,
		 * each edge is traversed twice (once forward, once reverse)
		 */
		if (isForwardEdge(e, g))
			outputEdge(e, g);

		return BFS_SUCCESS;
	}

	/** Invoked when a vertex is visited for the first time */
	BFSVisitorResult discover_vertex(const VertexT& v, const GraphT&)
	{
		/* declare vertex (GraphViz) */
		m_out << '\t' << v.kmer().c_str();

		m_out << " [";

		DepthMapConstIt depthIt = m_depthMap.find(v);
		assert(depthIt != m_depthMap.end());
		unsigned depth = depthIt->second;

	    m_out << "depth=" << depth;

		for (KmerPropertiesIt it = m_kmerProperties.begin();
			it != m_kmerProperties.end(); ++it) {
			if (it->second.find(v) != it->second.end())
				m_out << "," << it->first;
		}
		m_out << "];\n";

		return BFS_SUCCESS;
	}

	/**
	 * Invoked when an edge is traversed. (Each edge
	 * in the graph is traversed exactly once.)
	 */
	BFSVisitorResult non_tree_edge(const EdgeT& e, const GraphT& g)
	{
		/*
		 * since we use this visitor class with an undirected BFS,
		 * each edge is traversed twice (once forward, once reverse)
		 */
		if (isForwardEdge(e, g))
			outputEdge(e, g);
		return BFS_SUCCESS;
	}

	BFSVisitorResult examine_vertex(const VertexT& v, const GraphT&)
	{
		m_currentVertex = v;
		return BFS_SUCCESS;
	}

	/** Return if the current edge is a forward edge */
	bool isForwardEdge(const EdgeT& e, const GraphT& g)
	{
		assert(source(e, g) == m_currentVertex
			|| target(e, g) == m_currentVertex);
		return source(e, g) == m_currentVertex;
	}

	/** Output and edge */
	void outputEdge(const EdgeT& e, const GraphT& g)
	{
		const VertexT& u = source(e, g);
		const VertexT& v = target(e, g);

		/* declare edge (GraphViz) */
		m_out << '\t' << u.kmer().c_str() << " -> "
			<< v.kmer().c_str() << ";\n";
	}

protected:

	VertexT m_currentVertex;
	DepthMap m_depthMap;

	const VertexT& m_root;
	const unsigned m_maxDepth;
	const KmerProperties& m_kmerProperties;
	std::ostream& m_out;
};


#endif
