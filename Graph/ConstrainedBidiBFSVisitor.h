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
	V m_start;
	V m_goal;

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

	/** max number of edges to traverse during search */
	unsigned m_maxCost;
	/** number of edges traversed so far */
	unsigned m_cost;

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
	bool m_maxCostExceeded;
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
		unsigned maxCost,
		size_t memLimit
		) :
			m_graph(graph),
			m_start(start),
			m_goal(goal),
			m_maxPaths(maxPaths),
			m_minPathLength(minPathLength),
			m_maxPathLength(maxPathLength),
			m_maxBranches(maxBranches),
			m_maxCost(maxCost),
			m_cost(0),
			m_memLimit(memLimit),
			m_memCheckCounter(0),
			m_exceededMemLimit(false),
			m_peakActiveBranches(0),
			m_tooManyBranches(false),
			m_maxCostExceeded(false),
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
	void examine_vertex(const V& v, const G&, Direction dir)
	{
		std::cout << "visiting vertex: " << v << " from dir: " << dir << "\n";
	}

	void examine_edge(const E& e, const G& g, Direction dir)
	{
		V u = source(e, g);
		V v = target(e, g);
		std::cout << "visiting edge: (" << u << "," << v
			<< ") from dir: " << dir << "\n";
	}
#endif

	BFSVisitorResult discover_vertex(const V&, const G&,
			Direction, unsigned numActiveBranches)
	{
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
		if (m_cost >= m_maxCost) {
			m_maxCostExceeded = true;
			return ABORT_SEARCH;
		}
		m_cost++;

		if (!updateTargetDepth(e, g, dir))
			return SKIP_ELEMENT;

		return recordEdgeTraversal(e, g, dir);
	}

	BFSVisitorResult non_tree_edge(const E& e, const G& g, Direction dir)
	{
		if (m_cost >= m_maxCost) {
			m_maxCostExceeded = true;
			return ABORT_SEARCH;
		}
		m_cost++;

		return recordEdgeTraversal(e, g, dir);
	}

	BFSVisitorResult common_edge(const E& e, const G& g, Direction dir)
	{
		if (m_cost >= m_maxCost) {
			m_maxCostExceeded = true;
			return ABORT_SEARCH;
		}
		m_cost++;

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
		else if (m_maxCostExceeded)
			return MAX_COST_EXCEEDED;
		else if (m_exceededMemLimit)
			return EXCEEDED_MEM_LIMIT;

		PathSearchResult result = buildPaths();

		if (result == FOUND_PATH)
			pathsFound = m_pathsFound;

		return result;
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

	unsigned getSearchCost()
	{
		return m_cost;
	}

	size_t approxMemUsage()
	{
		return
			m_traversalGraph[FORWARD].approxMemSize() +
			m_traversalGraph[REVERSE].approxMemSize() +
			approxMemSize(m_depthMap[FORWARD]) +
			approxMemSize(m_depthMap[REVERSE]);
	}

	void getTraversalGraph(HashGraph<V>& traversalGraph)
	{
		typedef typename HashGraph<V>::vertex_iterator vertex_iterator;
		typedef typename HashGraph<V>::adjacency_iterator adjacency_iterator;

		Direction dir[] = { FORWARD, REVERSE };
		for (unsigned i = 0; i < 2; i++) {
			HashGraph<V>& g = m_traversalGraph[dir[i]];
			vertex_iterator vi, vi_end;
			boost::tie(vi, vi_end) = vertices(g);
			for(; vi != vi_end; vi++) {
				adjacency_iterator ai, ai_end;
				boost::tie(ai, ai_end) = adjacent_vertices(*vi, g);
				for(; ai != ai_end; ai++) {
					add_edge(*ai, *vi, traversalGraph);
				}
			}
		}
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
		const size_t MEM_COUNTER_ROLLOVER = 1000;
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

	PathSearchResult buildPaths()
	{
		PathSearchResult overallResult = NO_PATH;

		// m_pathsFound will already contain one sol'n
		// in the special case where start_kmer == goal_kmer
		if (!m_pathsFound.empty())
			overallResult = FOUND_PATH;

		typename EdgeSet::const_iterator i = m_commonEdges.begin();
		for (; i != m_commonEdges.end(); i++) {
			PathSearchResult result = buildPaths(*i);
			if (result == FOUND_PATH) {
				overallResult = FOUND_PATH;
			}
			else if (result != FOUND_PATH && result != NO_PATH) {
				// we have encountered a failure case
				// (e.g. TOO_MANY_PATHS)
				overallResult = result;
				break;
			}
		}
		return overallResult;
	}

	PathSearchResult buildPaths(const E& common_edge)
	{
		if (m_cost > m_maxCost) {
			m_maxCostExceeded = true;
			return MAX_COST_EXCEEDED;
		}

		V u = source(common_edge, m_graph);
		V v = target(common_edge, m_graph);

		// find paths from common edge to start vertex (forward traversal)

		unsigned maxPathsToStart = m_maxPaths - m_pathsFound.size();

		PathSearchResult resultCode;
	
		AllPathsSearchResult<V> leftResult = allPathsSearch(
			m_traversalGraph[FORWARD], u, m_start, maxPathsToStart,
			0, m_maxDepth[FORWARD], m_maxCost - m_cost);
		m_cost += leftResult.cost;
		resultCode = leftResult.resultCode;

		if (resultCode == FOUND_PATH) {

			// find paths from common edge to goal vertex (reverse traversal)

			unsigned maxPathsToGoal =
				(m_maxPaths - m_pathsFound.size()) / leftResult.paths.size();

			AllPathsSearchResult<V> rightResult =
				allPathsSearch(m_traversalGraph[REVERSE], v, m_goal,
				maxPathsToGoal, 0, m_maxDepth[REVERSE], m_maxCost - m_cost);
			m_cost += rightResult.cost;
			resultCode = rightResult.resultCode;

			if (resultCode == FOUND_PATH)
				resultCode = buildPaths(leftResult.paths, rightResult.paths);

		} // result == FOUND_PATH (common edge => start)

		if (resultCode == MAX_COST_EXCEEDED)
			m_maxCostExceeded = true;
		else if (resultCode == TOO_MANY_PATHS)
			m_tooManyPaths = true;

		return resultCode;
	}

	PathSearchResult buildPaths(const PathList& pathsToStart, const PathList& pathsToGoal)
	{
		bool addedPath = false;
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
				addedPath = true;
			}
		}
		return (addedPath ? FOUND_PATH : NO_PATH);
	}

};

#endif /* CONSTRAINED_BFS_VISITOR_H */
