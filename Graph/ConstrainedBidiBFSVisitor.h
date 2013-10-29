#ifndef CONSTRAINED_BIDI_BFS_VISITOR_H
#define CONSTRAINED_BIDI_BFS_VISITOR_H

#include "Common/UnorderedMap.h"
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

template <typename G>
class ConstrainedBidiBFSVisitor : public BidirectionalBFSVisitor<G>
{
public:

	typedef typename boost::graph_traits<G>::vertex_descriptor V;
	typedef typename boost::graph_traits<G>::edge_descriptor E;
	typedef unsigned short depth_t;

private:

	typedef unordered_map<V, depth_t> DepthMap;

	HashGraph<V> m_forwardTraversalGraph;
	HashGraph<V> m_reverseTraversalGraph;
	DepthMap m_depthMap;
	const V& m_start;
	const V& m_goal;
	unsigned m_maxPaths;
	depth_t m_minDepth;
	depth_t m_maxDepth;
	bool m_bFoundGoal;
	bool m_tooManyPaths;
	depth_t m_maxDepthVisited;

	std::vector< Path<V> > m_pathsFound;

public:

	ConstrainedBidiBFSVisitor(
		const V& start,
		const V& goal,
		unsigned maxPaths,
		depth_t minDepth,
		depth_t maxDepth) :
			m_start(start),
			m_goal(goal),
			m_maxPaths(maxPaths),
			m_minDepth(minDepth),
			m_maxDepth(maxDepth),
			m_bFoundGoal(false),
			m_tooManyPaths(false),
			m_maxDepthVisited(0)
	{
		// edge case: start == goal
		if (start == goal) {
			Path<V> path;
			path.push_back(start);
			m_pathsFound.push_back(path);
		}
	}

#if 0
	// for debugging
	void examine_vertex(const V& v, const G& g, Direction dir)
	{
		std::cout << "visiting vertex: " << v << "\n";
	}
#endif

	void record_edge_traversal(const E& e, const G& g, Direction dir)
	{
		V u = source(e, g);
		V v = target(e, g);

#if 0
		// for debugging
		std::cout << "visiting edge: (" << u << ", " << v << ")\n";
#endif

		// record history of traversal, so that we can retrace
		// paths from a common edge to start/goal.

		if (dir == FORWARD)
			add_edge(v, u, m_forwardTraversalGraph);
		else
			add_edge(u, v, m_reverseTraversalGraph);
	}

	BFSVisitorResult tree_edge(const E& e, const G& g, Direction dir)
	{
		record_edge_traversal(e, g, dir);

		V u = source(e, g);
		V v = target(e, g);

		if (m_depthMap[u] == m_maxDepth)
			return SKIP_ELEMENT;

		m_depthMap[v] = m_depthMap[u] + 1;
		if (m_depthMap[v] > m_maxDepthVisited)
			m_maxDepthVisited = m_depthMap[v];

		return SUCCESS;
	}

	BFSVisitorResult non_tree_edge(const E& e, const G& g, Direction dir)
	{
		record_edge_traversal(e, g, dir);
		return SUCCESS;
	}

	BFSVisitorResult common_edge(const E& e, const G& g)
	{
		if (m_pathsFound.size() == m_maxPaths) {
			m_tooManyPaths = true;
			return ABORT_SEARCH;
		}

		V u = source(e, g);
		V v = target(e, g);

		// find paths from common edge to start vertex (forward traversal)

		unsigned maxPathsToStart = m_maxPaths - m_pathsFound.size();
		std::vector< Path<V> > pathsToStart;
		PathSearchResult result = allPathsSearch(
			m_forwardTraversalGraph, u, m_start, maxPathsToStart,
			m_minDepth, m_maxDepth, pathsToStart);

		if (result == TOO_MANY_PATHS) {
			m_tooManyPaths = true;
			return ABORT_SEARCH;
		} else if (result == FOUND_PATH) {

			// find paths from common edge to goal vertex (reverse traversal)

			unsigned maxPathsToGoal = (m_maxPaths - m_pathsFound.size()) / pathsToStart.size();
			std::vector< Path<V> > pathsToGoal;
			result = allPathsSearch(m_forwardTraversalGraph, v, m_goal,
				maxPathsToGoal, m_minDepth, m_maxDepth, pathsToGoal);

			if (result == TOO_MANY_PATHS) {
				m_tooManyPaths = true;
				return ABORT_SEARCH;
			} else if (result == FOUND_PATH) {

				// merge paths from start to common edge
				// and common edge to goal

				for (unsigned i = 0; i < pathsToStart.size(); i++) {
					for(unsigned j = 0; j < pathsToGoal.size(); j++) {
						m_pathsFound.push_back(pathsToStart[i]);
						reverse(m_pathsFound.back().begin(), m_pathsFound.back().end());
						m_pathsFound.back().insert(m_pathsFound.back().end(),
							pathsToGoal[j].begin(), pathsToGoal[j].end());
					}
				}
			}

		}

		return SUCCESS;
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

	PathSearchResult pathsToGoal(std::vector< Path<V> >& pathsFound)
	{
		if (m_tooManyPaths) {
			return TOO_MANY_PATHS;
		} else if (m_pathsFound.empty()) {
			return NO_PATH;
		} else {
			pathsFound = m_pathsFound;
			return FOUND_PATH;
		}
	}

	depth_t getMaxDepthVisited()
	{
		return m_maxDepthVisited;
	}

};


#endif /* CONSTRAINED_BFS_VISITOR_H */
