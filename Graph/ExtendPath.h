
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
	/* longest branch of Bloom filter false positives */
	unsigned fpTrim;
	/* maximum length after extension */
	unsigned maxLen;
	/*
	 * if true, multiple incoming branches > trimLen
	 * will cause a path extension to halt
	 */
	bool lookBehind;
	/*
	 * If false, ignore incoming branches for the starting vertex.
	 * This is useful when when we are intentionally starting our
	 * path extension from a branching point.
	 */
    bool lookBehindStartVertex;

	/* constructor */
	ExtendPathParams() : trimLen(0), fpTrim(0), maxLen(NO_LIMIT), lookBehind(true),
		lookBehindStartVertex(true) {}
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
	if (depth >= depthLimit)
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
 * Return true if the given edge represents the start of a "true branch".
 * Roughly speaking, a path is a true branch if it has length >= trim
 * or terminates in a branching node, where a branching node is (recursively)
 * defined to be a node with either >= 2 incoming true branches or >= outgoing
 * true branches.
 */
template <class Graph>
static inline bool trueBranch(
	const typename boost::graph_traits<Graph>::edge_descriptor& e,
	unsigned depth, Direction dir, const Graph& g, unsigned trim,
	unsigned fpTrim, unordered_set<typename boost::graph_traits<Graph>::vertex_descriptor>& visited)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;

	typename boost::graph_traits<Graph>::out_edge_iterator oei, oei_end;
	typename boost::graph_traits<Graph>::in_edge_iterator iei, iei_end;

	const V& u = (dir == FORWARD) ? source(e, g) : target(e, g);
	const V& v = (dir == FORWARD) ? target(e, g) : source(e, g);

	/* branches with bubbles/cycles are considered true branches */
	if (visited.find(v) != visited.end())
		return true;

	if (depth >= trim)
		return true;

	visited.insert(v);

	if (dir == FORWARD) {
		for (boost::tie(oei, oei_end) = out_edges(v, g);
			oei != oei_end; ++oei) {
			if (trueBranch(*oei, depth+1, FORWARD, g, trim, fpTrim, visited))
				return true;
		}
		if (depth >= fpTrim || lookAhead(v, FORWARD, fpTrim, g)) {
			for (boost::tie(iei, iei_end) = in_edges(v, g);
				iei != iei_end; ++iei) {
				if (source(*iei, g) == u)
					continue;
				if (trueBranch(*iei, 0, REVERSE, g, trim, fpTrim, visited))
					return true;
			}
		}
	} else {
		assert(dir == REVERSE);
		for (boost::tie(iei, iei_end) = in_edges(v, g);
			iei != iei_end; ++iei) {
			if (trueBranch(*iei, depth+1, REVERSE, g, trim, fpTrim, visited))
				return true;
		}
		if (depth >= fpTrim || lookAhead(v, REVERSE, fpTrim, g)) {
			for (boost::tie(oei, oei_end) = out_edges(v, g);
				 oei != oei_end; ++oei) {
				if (target(*oei, g) == u)
					continue;
				if (trueBranch(*oei, 0, FORWARD, g, trim, fpTrim, visited))
					return true;
			}
		}
	}

	visited.erase(v);

	return false;
}

/**
 * Return true if the given edge represents the start of a "true branch".
 * Roughly speaking, a path is a true branch if it has length >= trim
 * or terminates in a branching node, where a branching node is (recursively)
 * defined to be a node with either >= 2 incoming true branches or >= outgoing
 * true branches.
 */
template <class Graph>
static inline bool trueBranch(
	const typename boost::graph_traits<Graph>::edge_descriptor& e,
	Direction dir, const Graph& g, unsigned trim, unsigned fpTrim)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;
	unordered_set<V> visited;
	return trueBranch(e, 0, dir, g, trim, fpTrim, visited);
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
	Direction dir, const BidirectionalGraph& g, unsigned trim, unsigned fpTrim)
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
			if (trueBranch(*oei, dir, g, trim, fpTrim))
				branchRoots.push_back(v);
		}
	} else {
		assert(dir == REVERSE);
		for (boost::tie(iei, iei_end) = in_edges(u, g);
			iei != iei_end; ++iei) {
			const V& v = source(*iei, g);
			if (trueBranch(*iei, dir, g, trim, fpTrim)) {
				branchRoots.push_back(v);
			}
		}
	}

	return branchRoots;
}

/**
 * Return the unique predecessor/successor of a given vertex. In
 * cases where a predecessor/successor does not exist (i.e. a dead end)
 * or is not unique (i.e. a branching point), return an result code indicating
 * why a unique successor could not be returned.
 */
template <class Graph>
static inline std::pair<typename boost::graph_traits<Graph>::vertex_descriptor,
	SingleExtensionResult>
successor(const typename boost::graph_traits<Graph>::vertex_descriptor& u,
	Direction dir, const Graph& g, unsigned trim, unsigned fpTrim)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;
	typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIt;
	typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIt;

	InEdgeIt iei, iei_end;
    OutEdgeIt oei, oei_end;

	V v;
	for (unsigned i = 0; true; i = (i == 0) ? 1 : std::min(trim, 2*i))
	{
		unsigned trueBranches = 0;

		if (dir == FORWARD) {
			for (boost::tie(oei, oei_end) = out_edges(u, g); oei != oei_end; ++oei) {
				if (trueBranch(*oei, FORWARD, g, i, fpTrim)) {
					v = target(*oei, g);
					++trueBranches;
					if (trueBranches >= 2)
						break;
				}
			}
		} else {
			assert(dir == REVERSE);
			for (boost::tie(iei, iei_end) = in_edges(u, g); iei != iei_end; ++iei) {
				if (trueBranch(*iei, REVERSE, g, i, fpTrim)) {
					v = source(*iei, g);
					++trueBranches;
					if (trueBranches >= 2)
						break;
				}
			}
		}

		if (trueBranches == 0)
			return std::make_pair(v, SE_DEAD_END);
		else if (trueBranches == 1)
			return std::make_pair(v, SE_EXTENDED);
		else if (i == trim)
			return std::make_pair(v, SE_BRANCHING_POINT);
	}

}

/**
 * Return true if the given vertex has more than one possible
 * predecessor/successor in the graph.
 */
template <class Graph>
static inline bool
ambiguous(const typename boost::graph_traits<Graph>::vertex_descriptor& u,
	Direction dir, const Graph& g, unsigned trim, unsigned fpTrim)
{
	return successor(u, dir, g, trim, fpTrim).second == SE_BRANCHING_POINT;
}

/**
 * Return true if the given vertex has more than one possible
 * predecessor/successor in the graph.
 *
 * @param expected always include this vertex in the set of possible
 * predecssors/successors, even if it is not a true branch.
 */
template <class Graph>
static inline bool
ambiguous(const typename boost::graph_traits<Graph>::vertex_descriptor& u,
	const typename boost::graph_traits<Graph>::vertex_descriptor& expected,
	Direction dir, const Graph& g, unsigned trim, unsigned fpTrim)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;

	V v;
	SingleExtensionResult result;

	boost::tie(v, result) = successor(u, dir, g, trim, fpTrim);

	return result == SE_BRANCHING_POINT
		|| (result == SE_EXTENDED && v != expected);
}

/**
 * Extend path by a single vertex, if there is a unique predecessor/successor
 * in the direction of extension.
 */
template <class Graph>
static inline SingleExtensionResult
extendPath(Path<typename boost::graph_traits<Graph>::vertex_descriptor>& path,
	Direction dir, const Graph& g, unsigned trim, unsigned fpTrim,
	bool lookBehind)
{
	assert(!path.empty());

	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;

	V t, v;
	SingleExtensionResult result;

	const V& head = (dir == FORWARD) ? path.back() : path.front();

	boost::tie(v, result) = successor(head, dir, g, trim, fpTrim);
	if (result != SE_EXTENDED)
		return result;

	if (lookBehind) {

		Direction otherDir = (dir == FORWARD) ? REVERSE : FORWARD;
		boost::tie(t, result) = successor(head, otherDir, g, trim, fpTrim);

		if (result == SE_BRANCHING_POINT)
			return result;

		/*
		 * Tricky: If our path was seeded on a tip, we want to stop the
		 * extension when we reach a branching point. We can detect that
		 * we are on tip if the previous path vertex does not match
		 * the expected predecessor `t`.
		 */
		if (path.size() > 1) {
			if (result == SE_DEAD_END) {
				/* no predecessors or all predecessors were tips */
				return SE_BRANCHING_POINT;
			} else {
				/* check if we are on a tip */
				assert(result == SE_EXTENDED);
				const V& prev = (dir == FORWARD) ?
					*(path.rbegin() + 1) : *(path.begin() + 1);
				if (prev != t)
					return SE_BRANCHING_POINT;
			}
		}

	}

	if (dir == FORWARD)
		path.push_back(v);
	else
		path.push_front(v);

	return SE_EXTENDED;
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
inline static
std::pair<typename boost::graph_traits<Graph>::vertex_descriptor, bool>
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
	bool tie = false;
	/* note: had to initialize to prevent compiler warnings */
	V longestBranch = u;
	if (dir == FORWARD) {
		for (boost::tie(oei, oei_end) = out_edges(u, g);
			 oei != oei_end; ++oei) {
			degree++;
			const V& v = target(*oei, g);
			size_t d = depth(v, dir, g) + 1;
			if (d > maxDepth) {
				maxDepth = d;
				longestBranch = v;
				tie = false;
			} else if (d == maxDepth && v < longestBranch) {
				/*
				 * make an arbitrary choice among branches
				 * of equal length using the vertex comparison
				 * operator (operator<).
				 */
				longestBranch = v;
				tie = true;
			}
		}
	} else {
		assert(dir == REVERSE);
		for (boost::tie(iei, iei_end) = in_edges(u, g);
			 iei != iei_end; ++iei) {
			degree++;
			const V& v = source(*iei, g);
			size_t d = depth(v, dir, g) + 1;
			if (d > maxDepth) {
				maxDepth = d;
				longestBranch = v;
				tie = false;
			} else if (d == maxDepth && v < longestBranch) {
				/*
				 * make an arbitrary choice among branches
				 * of equal length using the vertex comparison
				 * operator (operator<).
				 */
				longestBranch = v;
				tie = true;
			}
		}
	}
	assert(degree > 0);
	return std::make_pair(longestBranch, tie);
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
	bool cycle = false;
	bool lookBehind = params.lookBehindStartVertex;

	assert(!path.empty());
	while (result == SE_EXTENDED && path.size() < params.maxLen)
	{
		result = extendPath(path, dir, g, params.trimLen, params.fpTrim,
			lookBehind);

		if (result == SE_EXTENDED) {
			const V& head = (dir == FORWARD) ? path.back() : path.front();
			bool inserted;
			boost::tie(boost::tuples::ignore, inserted) = visited.insert(head);
			if (!inserted) {
				cycle = true;
				if (dir == FORWARD)
					path.pop_back();
				else
					path.pop_front();
				break;
			}
		}

		/* override `lookBehindStartVertex` after first extension */
		lookBehind = params.lookBehind;
	}

	if (path.size() > origPathLen) {
		if (cycle) {
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
		if (cycle) {
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
