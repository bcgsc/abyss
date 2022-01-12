
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
enum PathExtensionResultCode {
	/** stopped path extension at a vertex with multiple incoming branches */
	ER_AMBI_IN,
	/** stopped path extension at a vertex with multiple outgoing branches */
    ER_AMBI_OUT,
	/** stopped path extension at a vertex with no outgoing branches */
	ER_DEAD_END,
	/** stopped path extension after completing a cycle */
	ER_CYCLE,
	/** stopped path extension at caller-specified length limit */
	ER_LENGTH_LIMIT,
};

/**
 * Translate path extension result code to a string.
 */
static inline const char* pathExtensionResultStr(PathExtensionResultCode result)
{
	switch(result) {
	case ER_AMBI_IN:
		return "AMBI_IN";
	case ER_AMBI_OUT:
		return "AMBI_OUT";
	case ER_DEAD_END:
		return "DEAD_END";
	case ER_CYCLE:
		return "CYCLE";
	case ER_LENGTH_LIMIT:
		return "LENGTH_LIMIT";
	default:
		assert(false);
	}
	return "";
}

/** length of path extension (in vertices) and reason for stopping */
typedef std::pair<unsigned, PathExtensionResultCode> PathExtensionResult;

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
 * Return true if the given edge represents the beginning of a "true branch".
 *
 * A path is a true branch if it has length >= `trim` or terminates in a
 * branching node, where a branching node is (recursively) defined to be
 * a node with either >= 2 incoming true branches or >= 2 outgoing true branches.
 *
 * This method is similar to `lookAhead`, but it additionally changes traversal
 * direction when a dead-end is encountered.
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
		/*
		 * Note: The test for depth/lookAhead >= fpTrim before changing
		 * traversal direction is needed to deal with an X-shaped
		 * graph pattern that is frequently created by Bloom false positives.
		 * See the test for `trueBranch` in `ExtendPathTest.h` for an example.
		 */
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
		/*
		 * Note: The test for depth/lookAhead >= fpTrim before changing
		 * traversal direction is needed to deal with an X-shaped
		 * graph pattern that is frequently created by Bloom false positives.
		 * See the test for `trueBranch` in `ExtendPathTest.h` for an example.
		 */
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
	PathExtensionResultCode>
successor(const typename boost::graph_traits<Graph>::vertex_descriptor& u,
	Direction dir, const Graph& g, unsigned trim, unsigned fpTrim)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;
	typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIt;
	typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIt;

	InEdgeIt iei, iei_end;
    OutEdgeIt oei, oei_end;

	/* assign u to suppress uninitialized warning */
	V v = u;
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
			return std::make_pair(v, ER_DEAD_END);
		else if (trueBranches == 1)
			return std::make_pair(v, ER_LENGTH_LIMIT);
		else if (i == trim)
			return std::make_pair(v, ER_AMBI_OUT);
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
	return successor(u, dir, g, trim, fpTrim).second == ER_AMBI_OUT;
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
	PathExtensionResultCode result;

	boost::tie(v, result) = successor(u, dir, g, trim, fpTrim);

	return result == ER_AMBI_OUT || (result == ER_LENGTH_LIMIT && v != expected);
}

/**
 * Extend path by a single vertex, if there is a unique predecessor/successor
 * in the direction of extension.
 */
template <class Graph>
static inline PathExtensionResultCode
extendPathBySingleVertex(
	Path<typename boost::graph_traits<Graph>::vertex_descriptor>& path,
	Direction dir, const Graph& g, unsigned trim, unsigned fpTrim,
	bool lookBehind)
{
	assert(!path.empty());

	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;

	V t, v;
	PathExtensionResultCode result;

	const V& head = (dir == FORWARD) ? path.back() : path.front();

	if (lookBehind) {

		Direction otherDir = (dir == FORWARD) ? REVERSE : FORWARD;
		boost::tie(t, result) = successor(head, otherDir, g, trim, fpTrim);

		if (result == ER_AMBI_OUT)
			return ER_AMBI_IN;

		/*
		 * Tricky: If our path was seeded on a tip, we want to stop the
		 * extension when we reconnect to the graph. We can detect that
		 * we are on tip if we reach a branching point where the predecessor
		 * vertex in the path does not match the expected predecessor `t`.
		 */
		if (path.size() > 1) {
			if (result == ER_DEAD_END) {
				/* no predecessors or all predecessors were tips */
				return ER_AMBI_IN;
			} else {
				/* check if we are on a tip */
				assert(result == ER_LENGTH_LIMIT);
				const V& prev = (dir == FORWARD) ?
					*(path.rbegin() + 1) : *(path.begin() + 1);
				if (prev != t)
					return ER_AMBI_IN;
			}
		}

	}

	boost::tie(v, result) = successor(head, dir, g, trim, fpTrim);
	if (result != ER_LENGTH_LIMIT)
		return result;

	if (dir == FORWARD)
		path.push_back(v);
	else
		path.push_front(v);

	return ER_LENGTH_LIMIT;
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
	PathExtensionResultCode result = ER_DEAD_END;
	bool lookBehind = params.lookBehindStartVertex;

	assert(!path.empty());
	while (path.size() < params.maxLen)
	{
		result = extendPathBySingleVertex(path, dir, g,
			params.trimLen, params.fpTrim, lookBehind);

		if (result != ER_LENGTH_LIMIT)
			break;

		const V& head = (dir == FORWARD) ? path.back() : path.front();
		bool inserted;
		boost::tie(boost::tuples::ignore, inserted) = visited.insert(head);
		if (!inserted) {
			result = ER_CYCLE;
			if (dir == FORWARD)
				path.pop_back();
			else
				path.pop_front();
			break;
		}

		/* override `lookBehindStartVertex` after first extension */
		lookBehind = params.lookBehind;
	}

	if (params.maxLen != NO_LIMIT && path.size() == params.maxLen)
		result = ER_LENGTH_LIMIT;

	assert(path.size() >= origPathLen);
	unsigned extension = path.size() - origPathLen;

	/*
	 * Sanity check: If no length limit was imposed, we must have stopped
	 * the extension for some other reason (e.g. dead end)
	 */
	assert(params.maxLen != NO_LIMIT || result != ER_LENGTH_LIMIT);

	return std::make_pair(extension, result);
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
