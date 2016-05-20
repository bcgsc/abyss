
#ifndef _EXTENDPATH_H_
#define _EXTENDPATH_H_

#include "Graph/Path.h"
#include "Common/UnorderedSet.h"
#include "Common/UnorderedMap.h"
#include "Common/Hash.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <algorithm>

/**
 * Parameters for path extension.
 */
struct ExtendPathParams
{
	/* ignore branches shorter than or equal to this length */
	unsigned trimLen;
	/* maximum length after extension */
	unsigned maxLen;
	/*
	 * if true, multiple incoming branches > trimLen
	 * will cause a path extension to halt
	 */
	bool lookBehind;

	/* constructor */
	ExtendPathParams() : trimLen(0), maxLen(NO_LIMIT), lookBehind(true) {}
};

/**
 * The result of attempting to extend a path.
 */
enum PathExtensionResult {
	/** path could not be extended because of a dead end */
	DEAD_END,
	/** path could not be extended because of a branching point */
	BRANCHING_POINT,
	/** path could not be extended because of a cycle */
	CYCLE,
	/** path could not be extended because of caller-specified length limit */
	LENGTH_LIMIT,
	/** path was extended up to a dead end */
	EXTENDED_TO_DEAD_END,
	/** path was extended up to a branching point */
	EXTENDED_TO_BRANCHING_POINT,
	/** path was extended up to a cycle */
	EXTENDED_TO_CYCLE,
	/** path was extended up to caller-specified length limit */
	EXTENDED_TO_LENGTH_LIMIT
};

/**
 * Translate path extension result code to a string.
 */
static inline const char* pathExtensionResultStr(PathExtensionResult result)
{
	switch(result) {
	case DEAD_END:
		return "DEAD_END";
	case BRANCHING_POINT:
		return "BRANCHING_POINT";
	case CYCLE:
		return "CYCLE";
	case LENGTH_LIMIT:
		return "LENGTH_LIMIT";
	case EXTENDED_TO_DEAD_END:
		return "EXTENDED_TO_DEAD_END";
	case EXTENDED_TO_BRANCHING_POINT:
		return "EXTENDED_TO_BRANCHING_POINT";
	case EXTENDED_TO_CYCLE:
		return "EXTENDED_TO_CYCLE";
	case EXTENDED_TO_LENGTH_LIMIT:
		return "EXTENDED_TO_LENGTH_LIMIT";
	default:
		assert(false);
	}
}

/**
 * Return true if the path extension result code indicates
 * that the path was successfully extended by one or more nodes.
 */
static inline bool pathExtended(PathExtensionResult result)
{
	switch(result) {
	case DEAD_END:
	case BRANCHING_POINT:
	case CYCLE:
	case LENGTH_LIMIT:
		return false;
	default:
		return true;
	}
	assert(false);
}

/**
 * The result of attempting to extend a path
 * by a single neighbouring vertex.
 */
enum SingleExtensionResult {
	SE_DEAD_END,
	SE_BRANCHING_POINT,
	SE_EXTENDED
};

/**
 * Return true if there is a path of at least depthLimit vertices
 * that extends from given vertex u, otherwise return false.
 * Implemented using a bounded depth first search.
 *
 * @param start starting vertex for traversal
 * @param dir direction for traversal (FORWARD or REVERSE)
 * @param depth depth of current vertex u
 * @param depthLimit maximum depth to probe
 * @param g graph to use for traversal
 * @param visited vertices that have already been visited by the DFS
 * @return true if at least one path with length >= len
 * extends from v in direction dir, false otherwise
 */
template <class Graph>
static inline bool lookAhead(
	const typename boost::graph_traits<Graph>::vertex_descriptor& u,
	Direction dir, unsigned depth, unsigned depthLimit,
	unordered_set< typename boost::graph_traits<Graph>::vertex_descriptor,
	hash<typename boost::graph_traits<Graph>::vertex_descriptor> >& visited, const Graph& g)
{
    typedef typename boost::graph_traits<Graph>::vertex_descriptor V;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
    typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;

	OutEdgeIter oei, oei_end;
	InEdgeIter iei, iei_end;

	visited.insert(u);
	if (depth == depthLimit)
		return true;

	if (dir == FORWARD) {
		for (boost::tie(oei, oei_end) = out_edges(u, g);
			oei != oei_end; ++oei) {
			const V& v = target(*oei, g);
			if (visited.find(v) == visited.end()) {
				if(lookAhead(v, dir, depth+1, depthLimit, visited, g))
					return true;
			}
		}
	} else {
		assert(dir == REVERSE);
		for (boost::tie(iei, iei_end) = in_edges(u, g);
			 iei != iei_end; ++iei) {
			const V& v = source(*iei, g);
			if (visited.find(v) == visited.end()) {
				if(lookAhead(v, dir, depth+1, depthLimit, visited, g))
					return true;
			}
		}
	}

	return false;
}

/**
 * Return true if there is a path of at least 'depth' vertices
 * that extends from given vertex v, otherwise return false.
 * Implemented using a bounded depth first search.
 *
 * @param start starting vertex for traversal
 * @param dir direction for traversal (FORWARD or REVERSE)
 * @param depth length of path to test for
 * @param g graph to use for traversal
 * @return true if at least one path with length >= len
 * extends from v in direction dir, false otherwise
 */
template <class Graph>
static inline bool lookAhead(
	const typename boost::graph_traits<Graph>::vertex_descriptor& start,
	Direction dir, unsigned depth, const Graph& g)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;
	unordered_set< V, hash<V> > visited;
	return lookAhead(start, dir, 0, depth, visited, g);
}

/**
 * Return neighbour vertices that begin branches that are longer than trimLen.
 *
 * @param u root vertex
 * @param dir direction for neighbours (FORWARD or REVERSE)
 * @param g graph
 * @param trimLen ignore all branches less than or equal to this length
 * @return std::vector of neighbour vertices that start branches that are
 * greater than trimLen vertices in length
 */
template <class BidirectionalGraph>
static inline std::vector<typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor>
trueBranches(const typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor& u,
	Direction dir, const BidirectionalGraph& g, unsigned trimLen=0)
{
	typedef BidirectionalGraph G;
	typedef boost::graph_traits<G> graph_traits;
	typedef typename graph_traits::vertex_descriptor V;

	typename graph_traits::out_edge_iterator oei, oei_end;
	typename graph_traits::in_edge_iterator iei, iei_end;

	std::vector<V> branchRoots;

	if (dir == FORWARD) {
		for (boost::tie(oei, oei_end) = out_edges(u, g);
			oei != oei_end; ++oei) {
			const V& v = target(*oei, g);
			if (lookAhead(v, dir, trimLen, g))
				branchRoots.push_back(v);
		}
	} else {
		assert(dir == REVERSE);
		for (boost::tie(iei, iei_end) = in_edges(u, g);
			iei != iei_end; ++iei) {
			const V& v = source(*iei, g);
			if (lookAhead(v, dir, trimLen, g)) {
				branchRoots.push_back(v);
			}
		}
	}

	return branchRoots;
}

/**
 * Return the in/out degree of a vertex, disregarding branches
 * <= trimLen.
 *
 * @param u the vertex of interest
 * @param dir FORWARD for out degree, REVERSE for in degree
 * @param g the graph
 * @param trimLen branches less then or equal to this length
 * are ignored (unless they are the only option)
 * @return the in/out degree of u, ignoring branches <= trimLen
 */
template <typename Graph>
static inline unsigned trueDegree(
	const typename boost::graph_traits<Graph>::vertex_descriptor& u,
	Direction dir, const Graph& g, unsigned trimLen=0)
{
	typedef boost::graph_traits<Graph> graph_traits;
	typedef typename graph_traits::vertex_descriptor V;

	unsigned degree = (dir == FORWARD) ? out_degree(u, g) : in_degree(u, g);
	if (degree <= 1)
		return degree;

	std::vector<V> branches = trueBranches(u, dir, g, trimLen);
	/*
	 * Note: If branches.size() == 0, we know from above that
	 * we must have 2 or more short branches. This situation typically occurs
	 * near coverage gaps, where one of the branches is the correct choice.
	 * (During path extension, our heuristic is to choose the longest branch
	 * and to continue extending.)
	 */
	if (branches.size() == 0)
		return 1;

	return branches.size();
}

/**
 * Return the depth of the graph from the given source vertex,
 * i.e. the distance of the furthest node.  The depth is measured
 * by means of an exhaustive breadth first search.
 *
 * @param root starting vertex for traversal
 * @param dir direction for traversal (FORWARD or REVERSE)
 * @param g graph to use for traversal
 * @return the distance of the furthest vertex from root
 */
template <typename Graph>
static inline size_t depth(
	typename boost::graph_traits<Graph>::vertex_descriptor root,
	Direction dir, const Graph& g)
{
    typedef typename boost::graph_traits<Graph>::vertex_descriptor V;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
    typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;

	OutEdgeIter oei, oei_end;
	InEdgeIter iei, iei_end;

	unordered_set<V, hash<V> > visited;
	typedef unordered_map<V, size_t> DepthMap;
	DepthMap depthMap;
	std::deque<V> q;

	q.push_back(root);

	visited.insert(root);
	std::pair<typename DepthMap::iterator, bool> inserted =
		depthMap.insert(std::make_pair(root, 0));
	assert(inserted.second);

	size_t maxDepth = 0;
	while (!q.empty()) {
		V& u = q.front();
		visited.insert(u);
		typename DepthMap::const_iterator it = depthMap.find(u);
		assert(it != depthMap.end());
		size_t depth = it->second;
		if (depth > maxDepth)
			maxDepth = depth;
		if (dir == FORWARD) {
			for (boost::tie(oei, oei_end) = out_edges(u, g);
				oei != oei_end; ++oei) {
				V v = target(*oei, g);
				if (visited.find(v) == visited.end()) {
					visited.insert(v);
					std::pair<typename DepthMap::iterator, bool> inserted =
						depthMap.insert(std::make_pair(v, depth+1));
					assert(inserted.second);
					q.push_back(v);
				}
			}
		} else {
			assert(dir == REVERSE);
			for (boost::tie(iei, iei_end) = in_edges(u, g);
				iei != iei_end; ++iei) {
				V v = source(*iei, g);
				if (visited.find(v) == visited.end()) {
					visited.insert(v);
					std::pair<typename DepthMap::iterator, bool> inserted =
						depthMap.insert(std::make_pair(v, depth+1));
					assert(inserted.second);
					q.push_back(v);
				}
			}
		}
		q.pop_front();
	}

	return maxDepth;
}

/**
 * Return the neighbor vertex corresponding to the longest branch.  If there
 * are no neighbour vertices, an assertion will be thrown. If there
 * is a tie between branch lengths, the "winning" branch is chosen arbitrarily.
 *
 * @param u root vertex
 * @param dir direction of branches to consider (FORWARD or REVERSE)
 * @param g the graph
 * @return the vertex at the head of the longest branch
 */
template <typename Graph>
inline static typename boost::graph_traits<Graph>::vertex_descriptor
longestBranch(const typename boost::graph_traits<Graph>::vertex_descriptor& u,
	Direction dir, const Graph& g)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
    typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;

	OutEdgeIter oei, oei_end;
	InEdgeIter iei, iei_end;
	size_t maxDepth = 0;
	unsigned degree = 0;
	/* note: had to initialize to prevent compiler warnings */
	V longestBranch = u;
	if (dir == FORWARD) {
		for (boost::tie(oei, oei_end) = out_edges(u, g);
			 oei != oei_end; ++oei) {
			degree++;
			const V& v = target(*oei, g);
			size_t d = depth(v, dir, g);
			if (d >= maxDepth) {
				maxDepth = d;
				longestBranch = v;
			}
		}
	} else {
		assert(dir == REVERSE);
		for (boost::tie(iei, iei_end) = in_edges(u, g);
			 iei != iei_end; ++iei) {
			degree++;
			const V& v = source(*iei, g);
			size_t d = depth(v, dir, g);
			if (d >= maxDepth) {
				maxDepth = d;
				longestBranch = v;
			}
		}
	}
	assert(degree > 0);
	return longestBranch;
}

/**
 * If the given path has only one possible next/prev vertex in the graph,
 * append/prepend that vertex to the path.
 *
 * @param path the path to extend (a list of vertices)
 * @param dir direction of extension (FORWARD or REVERSE)
 * @param g the graph to use for traversal
 * @param params parameters controlling extension (e.g. trimLen)
 * @return PathExtensionResult: NO_EXTENSION, HIT_BRANCHING_POINT, or EXTENDED
 */
template <class BidirectionalGraph>
static inline SingleExtensionResult extendPathBySingleVertex(
	Path<typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor>& path,
	Direction dir, const BidirectionalGraph& g, const ExtendPathParams& params)
{
	typedef BidirectionalGraph G;
	typedef boost::graph_traits<G> graph_traits;
	typedef typename graph_traits::vertex_descriptor V;

	typename graph_traits::out_edge_iterator oei, oei_end;
	typename graph_traits::in_edge_iterator iei, iei_end;

	assert(dir == FORWARD || dir == REVERSE);

	V& u = (dir == FORWARD) ? path.back() : path.front();

	unsigned outDegree = (dir == FORWARD) ? out_degree(u, g) : in_degree(u, g);
	if (outDegree == 0) {
		return SE_DEAD_END;
	}

	unsigned inDegree = 0;
	if (params.lookBehind)
		inDegree = (dir == FORWARD) ? in_degree(u, g) : out_degree(u, g);

	if ((!params.lookBehind || inDegree <= 1) && outDegree == 1) {
		if (dir == FORWARD) {
			const V& v = target(*(out_edges(u, g).first), g);
			path.push_back(v);
		} else {
			assert(dir == REVERSE);
			const V& v = source(*(in_edges(u, g).first), g);
			path.push_front(v);
		}
		return SE_EXTENDED;
	}

	Direction otherDir = (dir == FORWARD) ? REVERSE : FORWARD;
	std::vector<V> longBranchesOut = trueBranches(u, dir, g, params.trimLen);
	std::vector<V> longBranchesIn;

	if (params.lookBehind) {
		longBranchesIn = trueBranches(u, otherDir, g, params.trimLen);
		/*
		 * Tricky: Make sure the path we are extending
		 * is treated as a valid incoming branch, even if it is less
		 * than trimLen. This can happen if we seeded the path on
		 * an error branch or a branch that has a coverage gap.
		 */
		if (path.size() > 1) {
			const V& predecessor = (dir == FORWARD) ?
				*(path.rbegin() + 1) : *(path.begin() + 1);
			if (std::find(longBranchesIn.begin(), longBranchesIn.end(),
				predecessor) == longBranchesIn.end()) {
				longBranchesIn.push_back(predecessor);
			}
		}
	}

	if ((params.lookBehind && longBranchesIn.size() > 1) ||
		longBranchesOut.size() > 1)
		return SE_BRANCHING_POINT;

	if (longBranchesOut.size() == 0) {
		/*
		 * If we have multiple branches that are shorter
		 * than the trim length then choose the longest one.
		 * (This type of situation usually occurs near
		 * coverage gaps.)
		 */
		V v = longestBranch(u, dir, g);
		if (dir == FORWARD)
			path.push_back(v);
		else
			path.push_front(v);

		return SE_EXTENDED;
	}

	if (dir == FORWARD)
		path.push_back(longBranchesOut.front());
	else
		path.push_front(longBranchesOut.front());

	return SE_EXTENDED;
}

/**
 * Extend a path up to the next branching point in the graph.
 *
 * @param path path to extend (modified by this function)
 * @param dir direction to extend path (FORWARD or REVERSE)
 * @param g graph in which to perform the extension
 * @param visited set of previously visited vertices (used
 * to detect cycles in the de Bruijn graph)
 * @param params parameters controlling extension (e.g. trimLen)
 * @return PathExtensionResult: NO_EXTENSION, HIT_BRANCHING_POINT,
 * or EXTENDED.
 */
template <class BidirectionalGraph>
static inline PathExtensionResult extendPath(
	Path<typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor>& path,
	Direction dir, const BidirectionalGraph& g,
	unordered_set<typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor>& visited,
	const ExtendPathParams& params)
{
	typedef BidirectionalGraph G;
	typedef boost::graph_traits<G> graph_traits;
	typedef typename graph_traits::vertex_descriptor V;
	typename graph_traits::out_edge_iterator oei, oei_end;
	typename graph_traits::in_edge_iterator iei, iei_end;

	assert(path.size() > 0);
	size_t origPathLen = path.size();

	if (path.size() != NO_LIMIT && path.size() >= params.maxLen)
		return LENGTH_LIMIT;

	SingleExtensionResult result = SE_EXTENDED;
	bool detectedCycle = false;

	while (result == SE_EXTENDED && !detectedCycle &&
		path.size() < params.maxLen)
	{
		result = extendPathBySingleVertex(path, dir, g, params);
		if (result == SE_EXTENDED) {
			std::pair<typename unordered_set<V>::iterator,bool> inserted;
			if (dir == FORWARD) {
				inserted = visited.insert(path.back());
			} else {
				assert(dir == REVERSE);
				inserted = visited.insert(path.front());
			}
			if (!inserted.second)
				detectedCycle = true;
		}
	}

	/** the last kmer we added is a repeat, so remove it */
	if (detectedCycle) {
		if (dir == FORWARD) {
			path.pop_back();
		} else {
			assert(dir == REVERSE);
			path.pop_front();
		}
	}

	if (path.size() > origPathLen) {
		if (detectedCycle) {
			return EXTENDED_TO_CYCLE;
		} else if (result == SE_DEAD_END) {
			return EXTENDED_TO_DEAD_END;
		} else if (result == SE_BRANCHING_POINT) {
			return EXTENDED_TO_BRANCHING_POINT;
		} else {
			assert(result == SE_EXTENDED &&
				path.size() == params.maxLen);
			return EXTENDED_TO_LENGTH_LIMIT;
		}
	} else {
		assert(path.size() == origPathLen);
		if (detectedCycle) {
			return CYCLE;
		} else if (result == SE_DEAD_END) {
			return DEAD_END;
		} else if (result == SE_BRANCHING_POINT) {
			return BRANCHING_POINT;
		} else {
			assert(origPathLen >= params.maxLen);
			return LENGTH_LIMIT;
		}
	}
}

/**
 * Extend a path up to the next branching point in the graph.
 *
 * @param path path to extend (modified by this function)
 * @param dir direction to extend path (FORWARD or REVERSE)
 * @param g graph in which to perform the extension
 * @param params parameters controlling extension (e.g. trimLen)
 * @return PathExtensionResult: NO_EXTENSION, HIT_BRANCHING_POINT,
 * or EXTENDED.
 */
template <class BidirectionalGraph>
PathExtensionResult extendPath(
	Path<typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor>& path,
	Direction dir, const BidirectionalGraph& g, const ExtendPathParams& params)
{
	typedef typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor V;

	/* track visited nodes to avoid infinite traversal of cycles */
	unordered_set<V> visited;
	visited.insert(path.begin(), path.end());

	return extendPath(path, dir, g, visited, params);
}

/**
 * Extend a path up to the next branching point in the graph.
 *
 * @param path path to extend (modified by this function)
 * @param dir direction to extend path (FORWARD or REVERSE)
 * @param g graph in which to perform the extension
 * @return PathExtensionResult: NO_EXTENSION, HIT_BRANCHING_POINT,
 * or EXTENDED.
 */
template <class BidirectionalGraph>
PathExtensionResult extendPath(
	Path<typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor>& path,
	Direction dir, const BidirectionalGraph& g)
{
	/* default extension params */
	ExtendPathParams params;
	return extendPath(path, dir, g, params);
}

#endif
