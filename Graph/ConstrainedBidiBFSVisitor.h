#ifndef CONSTRAINED_BIDI_BFS_VISITOR_H
#define CONSTRAINED_BIDI_BFS_VISITOR_H

#include "Common/UnorderedMap.h"
#include "Common/UnorderedSet.h"
#include "Common/IOUtil.h"
#include "Graph/Path.h"
#include "Graph/HashGraph.h"
#include "Graph/BidirectionalBFSVisitor.h"
#include "Graph/AllPathsSearch.h"
#include <boost/graph/graph_traits.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <boost/functional/hash.hpp> // for boost::hash_combine

template <typename G>
class ConstrainedBidiBFSVisitor : public BidirectionalBFSVisitor<G>
{

protected:

	typedef typename boost::graph_traits<G>::vertex_descriptor V;
	typedef typename boost::graph_traits<G>::edge_descriptor E;
	typedef unsigned short depth_t;
	typedef std::vector< Path<V> > PathList;
	typedef unordered_map<V, depth_t> DepthMap;

	struct EdgeHash {
		std::size_t operator()(const E& e) const {
			std::size_t hash = 0;
			boost::hash_combine(hash, std::min(e.m_source, e.m_target));
			boost::hash_combine(hash, std::max(e.m_source, e.m_target));
			return hash;
		}
	};

	typedef unordered_set<E, EdgeHash> EdgeSet;

	const G& m_graph;
	const V& m_start;
	const V& m_goal;
	unsigned m_maxPaths;

	/** records history of forward/reverse traversals */
	HashGraph<V> m_traversalGraph[2];

	/** records depth of vertices during forward/reverse traversal */
	DepthMap m_depthMap[2];

	/** depth limits for forward/reverse traversal */
	depth_t m_maxDepth[2];

	/** edges that connect the forward and reverse traversals */
	EdgeSet m_commonEdges;

	depth_t m_minPathLength;
	depth_t m_maxPathLength;
	bool m_tooManyPaths;

	PathList m_pathsFound;

public:

	ConstrainedBidiBFSVisitor(
		const G& graph,
		const V& start,
		const V& goal,
		unsigned maxPaths,
		depth_t minPathLength,
		depth_t maxPathLength) :
			m_graph(graph),
			m_start(start),
			m_goal(goal),
			m_maxPaths(maxPaths),
			m_minPathLength(minPathLength),
			m_maxPathLength(maxPathLength),
			m_tooManyPaths(false)
	{

		depth_t maxDepth = maxPathLength - 1;
		m_maxDepth[FORWARD] = maxDepth / 2 + maxDepth % 2;
		m_maxDepth[REVERSE] = maxDepth / 2;

		// special case
		if (start == goal) {
			Path<V> path;
			path.push_back(start);
			m_pathsFound.push_back(path);
		}
	}

#if 1
	// for debugging
	void examine_vertex(const V& v, const G& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(g);
		std::cout << "visiting vertex: " << v << " from dir: " << dir << "\n";
	}

	void examine_edge(const E& e, const G& g, Direction dir)
	{
		SUPPRESS_UNUSED_WARNING(g);
		std::cout << "visiting edge: " << e << " from dir: " << dir << "\n";
	}
#endif

	BFSVisitorResult tree_edge(const E& e, const G& g, Direction dir)
	{
		recordEdgeTraversal(e, g, dir);

		V u = source(e, g);
		V v = target(e, g);

		if (!updateVertexDepth(u, v, dir))
			return SKIP_ELEMENT;

		return SUCCESS;
	}

	BFSVisitorResult non_tree_edge(const E& e, const G& g, Direction dir)
	{
		recordEdgeTraversal(e, g, dir);
		return SUCCESS;
	}

	BFSVisitorResult common_tree_edge(const E& e, const G& g, Direction dir)
	{
		V u = source(e, g);
		V v = target(e, g);

		if (!updateVertexDepth(u, v, dir))
			return SKIP_ELEMENT;

		return recordCommonEdge(e);
	}

	BFSVisitorResult common_non_tree_edge(const E& e, const G& g)
	{
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

		buildPaths();

		if (m_tooManyPaths) {
			return TOO_MANY_PATHS;
		}
		else if (m_pathsFound.empty()) {
			return NO_PATH;
		} else {
			pathsFound = m_pathsFound;
			return FOUND_PATH;
		}
	}

protected:

	BFSVisitorResult recordCommonEdge(const E& e)
	{
		m_commonEdges.insert(e);
		if (m_commonEdges.size() > m_maxPaths) {
			m_tooManyPaths = true;
			return ABORT_SEARCH;
		}
		return SUCCESS;
	}

	/**
	 * Record history of edge traversal, so that we can retrace
	 * paths from a common edge to start/goal.
	 */
	void recordEdgeTraversal(const E& e, const G& g, Direction dir)
	{
		V u = source(e, g);
		V v = target(e, g);

		if (dir == FORWARD)
			add_edge(v, u, m_traversalGraph[FORWARD]);
		else
			add_edge(u, v, m_traversalGraph[REVERSE]);
	}

	/**
	 * Record the depth of a newly visited vertex.
	 * @return true if the vertex is visitable is less than the max
	 * depth limit false otherwise.
	 */
	bool updateVertexDepth(const V& u, const V& v, Direction dir)
	{
		const V* parent;
		const V* child;

		if (dir == FORWARD) {
			parent = &u;
			child = &v;
		} else {
			parent = &v;
			child = &u;
		}

		depth_t parentDepth = m_depthMap[dir][*parent];
		if (parentDepth == m_maxDepth[dir])
			return false;

		m_depthMap[dir][*child] = parentDepth + 1;
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
