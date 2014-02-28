#ifndef CONSTRAINED_BIDI_BFS_VISITOR_H
#define CONSTRAINED_BIDI_BFS_VISITOR_H

#include "Common/UnorderedMap.h"
#include "Common/UnorderedSet.h"
#include "Common/IOUtil.h"
#include "Graph/Path.h"
#include "Graph/HashGraph.h"
#include "Graph/BidirectionalBFSVisitor.h"
#include "Graph/AllPathsSearch.h"
#include "Common/MemUtils.h"
#include <boost/graph/graph_traits.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

template <typename G>
class ConstrainedBidiBFSVisitor : public BidirectionalBFSVisitor<G>
{

protected:

	static const float MEM_CHECK_FREQ = 0.001;
	static const size_t MEM_COUNTER_ROLLOVER = 1 / 0.001;

	typedef typename boost::graph_traits<G>::vertex_descriptor V;
	typedef typename boost::graph_traits<G>::edge_descriptor E;
	typedef unsigned short depth_t;
	typedef std::vector< Path<V> > PathList;
	typedef unordered_map<V, depth_t, hash<V> > DepthMap;

	struct EdgeHash {
		const G& m_g;
		EdgeHash(const G& g) : m_g(g) { }
		std::size_t operator()(const E& e) const {
			V u = source(e, m_g);
			V v = target(e, m_g);
			return hash<V>()(u) ^ hash<V>()(v);
		}
	};

	typedef unordered_set<E, EdgeHash> EdgeSet;

	const G& m_graph;
	const V& m_start;
	const V& m_goal;

	/** maximum number of paths to discover before aborting search */
	unsigned m_maxPaths;

	/** records history of forward/reverse traversals */
	HashGraph<V> m_traversalGraph[2];

	/** records depth of vertices during forward/reverse traversal */
	DepthMap m_depthMap[2];

	/** depth limits for forward/reverse traversal */
	depth_t m_maxDepth[2];

	/** max depth for forward/reverse traversal */
	depth_t m_maxDepthVisited[2];

	depth_t m_minPathLength;
	depth_t m_maxPathLength;

	/** maximum number of frontier nodes allowed at any given
	  * time during forward/reverse traversal */
	unsigned m_maxBranches;

	/** memory limit for graph search */
	size_t m_memLimit;
	/** controls frequency of memory limit checks */
	size_t m_memCheckCounter;
	/** true if we have exceeded the memory limit */
	bool m_exceededMemLimit;

	/** the max number of frontier nodes we had at any time
	  * during forward/reverse traversal (up to a limit
	  * of m_maxBranches) */
	unsigned m_peakActiveBranches;

	bool m_tooManyBranches;
	bool m_tooManyPaths;

	unsigned long long m_numNodesVisited;

	/** edges that connect the forward and reverse traversals */
	EdgeSet m_commonEdges;

	PathList m_pathsFound;

public:

	ConstrainedBidiBFSVisitor(
		const G& graph,
		const V& start,
		const V& goal,
		unsigned maxPaths,
		depth_t minPathLength,
		depth_t maxPathLength,
		unsigned maxBranches,
		size_t memLimit
		) :
			m_graph(graph),
			m_start(start),
			m_goal(goal),
			m_maxPaths(maxPaths),
			m_minPathLength(minPathLength),
			m_maxPathLength(maxPathLength),
			m_maxBranches(maxBranches),
			m_memLimit(memLimit),
			m_memCheckCounter(0),
			m_exceededMemLimit(false),
			m_peakActiveBranches(0),
			m_tooManyBranches(false),
			m_tooManyPaths(false),
			m_numNodesVisited(0),
			m_commonEdges(m_maxPaths, EdgeHash(m_graph))
	{

		depth_t maxDepth = maxPathLength - 1;
		m_maxDepth[FORWARD] = maxDepth / 2 + maxDepth % 2;
		m_maxDepth[REVERSE] = maxDepth / 2;

		m_maxDepthVisited[FORWARD] = 0;
		m_maxDepthVisited[REVERSE] = 0;

		// special case
		if (start == goal && 1 >= m_minPathLength) {
			Path<V> path;
			path.push_back(start);
			m_pathsFound.push_back(path);
		}
	}

#if 0
	// for debugging
	void examine_vertex(const V& v, const G& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(g);
		std::cout << "visiting vertex: " << v << " from dir: " << dir << "\n";
	}

	void examine_edge(const E& e, const G& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(g);
		V u = source(e, g);
		V v = target(e, g);
		std::cout << "visiting edge: (" << u << "," << v
			<< ") from dir: " << dir << "\n";
	}
#endif

	BFSVisitorResult discover_vertex(const V& v, const G& g,
			Direction dir, unsigned numActiveBranches)
	{
		SUPPRESS_UNUSED_WARNING(v);
		SUPPRESS_UNUSED_WARNING(g);
		SUPPRESS_UNUSED_WARNING(dir);

		if (m_maxBranches != NO_LIMIT &&
			numActiveBranches >= m_maxBranches) {
			m_tooManyBranches = true;
			return ABORT_SEARCH;
		}

		m_numNodesVisited++;

		// include new branch started by vertex v
		numActiveBranches++;
		if (numActiveBranches > m_peakActiveBranches)
			m_peakActiveBranches = numActiveBranches;

		return SUCCESS;
	}

	BFSVisitorResult tree_edge(const E& e, const G& g, Direction dir)
	{
		if (!updateTargetDepth(e, g, dir))
			return SKIP_ELEMENT;

		return recordEdgeTraversal(e, g, dir);
	}

	BFSVisitorResult non_tree_edge(const E& e, const G& g, Direction dir)
	{
		return recordEdgeTraversal(e, g, dir);
	}

	BFSVisitorResult common_edge(const E& e, const G& g, Direction dir)
	{
		V u = source(e, g);
		V v = target(e, g);

		const V& parent = (dir == FORWARD) ? u : v;

		if (m_depthMap[dir][parent] >= m_maxDepth[dir])
			return SKIP_ELEMENT;

		return recordCommonEdge(e);
	}

	PathSearchResult uniquePathToGoal(Path<V>& path)
	{
		std::vector< Path<V> > paths;
		PathSearchResult result = pathsToGoal(paths);
		if (paths.size() > 1) {
			return TOO_MANY_PATHS;
		} else if (result == FOUND_PATH && paths.size() == 1) {
			path = paths[0];
			return FOUND_PATH;
		}
		return result;
	}

	PathSearchResult pathsToGoal(PathList& pathsFound)
	{
		if (m_tooManyPaths)
			return TOO_MANY_PATHS;
		else if (m_tooManyBranches)
			return TOO_MANY_BRANCHES;
		else if (m_exceededMemLimit)
			return EXCEEDED_MEM_LIMIT;

		buildPaths();

		if (m_tooManyPaths) {
			return TOO_MANY_PATHS;
		} else if (m_pathsFound.empty()) {
			return NO_PATH;
		} else {
			pathsFound = m_pathsFound;
			return FOUND_PATH;
		}
	}

	depth_t getMaxDepthVisited(Direction dir)
	{
		return m_maxDepthVisited[dir];
	}

	unsigned getMaxActiveBranches()
	{
		return m_peakActiveBranches;
	}

	unsigned long long getNumNodesVisited()
	{
		return m_numNodesVisited;
	}

	size_t approxMemUsage()
	{
		return
			m_traversalGraph[FORWARD].approxMemSize() +
			m_traversalGraph[REVERSE].approxMemSize() +
			approxMemSize(m_depthMap[FORWARD]) +
			approxMemSize(m_depthMap[REVERSE]);
	}

protected:

	BFSVisitorResult recordCommonEdge(const E& e)
	{
		m_commonEdges.insert(e);
		if (m_maxPaths != NO_LIMIT &&
			m_commonEdges.size() > m_maxPaths) {
			m_tooManyPaths = true;
			return ABORT_SEARCH;
		}

		/**
		 * Tricky point:
		 *
		 * Recording the common edges in the both the
		 * forward and reverse traversal histories
		 * is necessary for edge cases where forward
		 * or reverse traversals are limited to
		 * a depth of zero. (In other words,
		 * the traversal graph has zero edges
		 * and exactly one vertex which is either
		 * the start or the goal vertex.)
		 *
		 * I cannot find a way to add a vertex
		 * with a specific vertex_descriptor
		 * to a graph using the Boost graph API.
		 * The only way seems to be creating an edge
		 * that has the given vertex_descriptor
		 * as the source or target.
		 */

		BFSVisitorResult result = recordEdgeTraversal(e, m_graph, FORWARD);
		if (result != SUCCESS)
			return result;

		return recordEdgeTraversal(e, m_graph, REVERSE);
	}

	BFSVisitorResult checkMemLimit()
	{
		m_memCheckCounter++;
		if (m_memCheckCounter >= MEM_COUNTER_ROLLOVER) {
			m_memCheckCounter = 0;
			if (approxMemUsage() > m_memLimit) {
				m_exceededMemLimit = true;
				return ABORT_SEARCH;
			}
		}
		return SUCCESS;
	}

	/**
	 * Record history of edge traversal, so that we can retrace
	 * paths from a common edge to start/goal.
	 */
	BFSVisitorResult recordEdgeTraversal(const E& e, const G& g, Direction dir)
	{
		BFSVisitorResult result = checkMemLimit();
		if (result != SUCCESS)
			return result;

		V u = source(e, g);
		V v = target(e, g);

		if (dir == FORWARD)
			add_edge(v, u, m_traversalGraph[FORWARD]);
		else
			add_edge(u, v, m_traversalGraph[REVERSE]);

		return result;
	}

	/**
	 * Record the depth of a newly visited vertex.
	 * @return true if the vertex is visitable is less than the max
	 * depth limit false otherwise.
	 */
	bool updateTargetDepth(const E& e, const G& g, Direction dir)
	{
		const V& parent = (dir == FORWARD) ? source(e, g) : target(e, g);
		const V& child = (dir == FORWARD) ? target(e, g) : source(e, g);

		depth_t parentDepth = m_depthMap[dir][parent];
		if (parentDepth == m_maxDepth[dir])
			return false;

		depth_t childDepth = parentDepth + 1;
		m_depthMap[dir][child] = childDepth;

		if (childDepth > m_maxDepthVisited[dir])
			m_maxDepthVisited[dir] = childDepth;

		return true;
	}

	void buildPaths()
	{
		typename EdgeSet::const_iterator i = m_commonEdges.begin();
		for (; i != m_commonEdges.end(); i++) {
			if (buildPaths(*i) == TOO_MANY_PATHS)
				break;
		}
	}

	PathSearchResult buildPaths(const E& common_edge)
	{
		V u = source(common_edge, m_graph);
		V v = target(common_edge, m_graph);

		// find paths from common edge to start vertex (forward traversal)

		unsigned maxPathsToStart = m_maxPaths - m_pathsFound.size();

		PathList pathsToStart;
		PathSearchResult result = allPathsSearch(
			m_traversalGraph[FORWARD], u, m_start, maxPathsToStart,
			0, m_maxDepth[FORWARD], pathsToStart);

		if (result == FOUND_PATH) {

			// find paths from common edge to goal vertex (reverse traversal)

			unsigned maxPathsToGoal =
				(m_maxPaths - m_pathsFound.size()) / pathsToStart.size();

			PathList pathsToGoal;
			result = allPathsSearch(m_traversalGraph[REVERSE], v, m_goal,
				maxPathsToGoal, 0, m_maxDepth[REVERSE], pathsToGoal);

			if (result == FOUND_PATH)
				buildPaths(pathsToStart, pathsToGoal);

		} // result == FOUND_PATH (common edge => start)

		if (result == TOO_MANY_PATHS)
			m_tooManyPaths = true;

		return result;
	}

	void buildPaths(const PathList& pathsToStart, const PathList& pathsToGoal)
	{
		typename PathList::const_iterator pathToStart = pathsToStart.begin();
		for (; pathToStart != pathsToStart.end(); pathToStart++) {
			typename PathList::const_iterator pathToGoal = pathsToGoal.begin();
			for(; pathToGoal != pathsToGoal.end(); pathToGoal++) {
				if (pathToStart->size() + pathToGoal->size() < m_minPathLength ||
					pathToStart->size() + pathToGoal->size() > m_maxPathLength)
					continue;
				m_pathsFound.push_back(*pathToStart);
				Path<V>& mergedPath = m_pathsFound.back();
				reverse(mergedPath.begin(), mergedPath.end());
				m_pathsFound.back().insert(mergedPath.end(),
					pathToGoal->begin(), pathToGoal->end());
			}
		}
	}

};

#endif /* CONSTRAINED_BFS_VISITOR_H */
