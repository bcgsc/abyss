#ifndef CONSTRAINED_BFS_VISITOR_H
#define CONSTRAINED_BFS_VISITOR_H

#include "Common/UnorderedMap.h"
#include "Common/IOUtil.h"
#include "Graph/BreadthFirstSearch.h"
#include "Graph/DefaultColorMap.h"
#include "Graph/Path.h"
#include "Graph/HashGraph.h"
#include "Graph/AllPathsSearch.h"
#include <boost/graph/graph_traits.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

template <typename G>
class ConstrainedBFSVisitor : public DefaultBFSVisitor<G>
{
public:

	typedef typename boost::graph_traits<G>::vertex_descriptor V;
	typedef typename boost::graph_traits<G>::edge_descriptor E;
	typedef unsigned short depth_t;

private:

	typedef std::vector<V> Predecessors;
	typedef unordered_map<V, unsigned, hash<V> > OutDegreeMap;
	typedef unordered_map<V, depth_t, hash<V> > DepthMap;

	HashGraph<V> m_traversalGraph;
	OutDegreeMap m_outDegree;
	DepthMap m_depthMap;
	V m_start;
	V m_goal;
	depth_t m_minDepth;
	depth_t m_maxDepth;
	unsigned m_maxBranches;
	DefaultColorMap<G>& m_colorMap;
	bool m_bFoundGoal;
	depth_t m_maxDepthVisited;
	unsigned m_branches;
	bool m_tooManyBranches;

public:

	ConstrainedBFSVisitor(
		const V& start,
		const V& goal,
		depth_t minDepth,
		depth_t maxDepth,
		unsigned maxBranches,
		DefaultColorMap<G>& colorMap) :
			m_start(start),
			m_goal(goal),
			m_minDepth(minDepth),
			m_maxDepth(maxDepth),
			m_maxBranches(maxBranches),
			m_colorMap(colorMap),
			m_bFoundGoal(false),
			m_maxDepthVisited(0),
			m_branches(1),
			m_tooManyBranches(false) {}

#if 0
	// for debugging
	BFSVisitorResult examine_vertex(const V& v, const G& g)
	{
		std::cout << "visiting vertex: " << v << "\n";
		return BFS_SUCCESS;
	}
#endif

	BFSVisitorResult examine_edge(const E& e, const G& g)
	{

		V u = source(e, g);
		V v = target(e, g);

#if 0
		// for debugging
		std::cout << "visiting edge: (" << u << ", " << v << ")\n";
#endif

		if (m_tooManyBranches) {
			put(m_colorMap, v, boost::black_color);
			return BFS_ABORT_SEARCH;
		}

		if (v == m_goal)
			m_bFoundGoal = true;

		// record history of traversal, so that we can trace
		// backwards from goal to start in pathsToGoal()

		add_edge(v, u, m_traversalGraph);

		// track depth of nodes and go no deeper than m_maxDepth

		if (get(m_colorMap, v) == boost::white_color) // tree edge
			m_depthMap[v] = m_depthMap[u] + 1;

		if (m_depthMap[v] >= m_maxDepth)
			put(m_colorMap, v, boost::black_color);

		if (m_depthMap[v] > m_maxDepthVisited)
			m_maxDepthVisited = m_depthMap[v];

		// track number of branches and abort if we exceed m_maxBranches

		if (m_maxBranches != NO_LIMIT && ++m_outDegree[u] > 1)
			m_branches++;

		if (m_maxBranches != NO_LIMIT && m_branches > m_maxBranches) {
			m_tooManyBranches = true;
			put(m_colorMap, v, boost::black_color);
		}

		return BFS_SUCCESS;
	}

	AllPathsSearchResult<V> uniquePathToGoal()
	{
		AllPathsSearchResult<V> result = pathsToGoal(1);
		return result;
	}

	AllPathsSearchResult<V> pathsToGoal(unsigned maxPaths)
	{
		AllPathsSearchResult<V> result;

		if (m_tooManyBranches) {
			result.resultCode = TOO_MANY_BRANCHES;
			return result;
		}
		else if (!m_bFoundGoal) {
			result.resultCode = NO_PATH;
			return result;
		}

		result = allPathsSearch(m_traversalGraph, m_goal, m_start,
			maxPaths, m_minDepth, m_maxDepth, NO_LIMIT);

		if (result.resultCode == FOUND_PATH) {
			for (unsigned i = 0; i < result.paths.size(); i++)
				reverse(result.paths[i].begin(), result.paths[i].end());
		}

		return result;
	}

	depth_t getMaxDepthVisited()
	{
		return m_maxDepthVisited;
	}

};


#endif /* CONSTRAINED_BFS_VISITOR_H */
