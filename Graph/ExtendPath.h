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

/**
 * The result of attempting to extend a path.
 */
enum PathExtensionResult {
	DEAD_END,
	BRANCHING_POINT,
	CYCLE,
	EXTENDED_TO_DEAD_END,
	EXTENDED_TO_BRANCHING_POINT,
	EXTENDED_TO_CYCLE
};

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
 * Return true if there is a path of at least 'depth' vertices
 * that extends from given vertex v, otherwise return false.
 * Implemented using a bounded breadth first search.
 *
 * @param start starting vertex for traversal
 * @param dir direction for traversal (FORWARD or REVERSE)
 * @param depth length of path to test for
 * @param g graph to use for traversal
 * @return true if at least one path with length >= len
 * extends from v in direction dir, false otherwise
 */
template <class BidirectionalGraph>
static inline bool lookAhead(
	typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor start,
	Direction dir, unsigned depth, const BidirectionalGraph& g)
{
    typedef typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor V;
    typedef typename boost::graph_traits<BidirectionalGraph>::out_edge_iterator OutEdgeIter;
    typedef typename boost::graph_traits<BidirectionalGraph>::in_edge_iterator InEdgeIter;

	OutEdgeIter oei, oei_end;
	InEdgeIter iei, iei_end;

	unordered_set<V, hash<V> > visited;
	typedef unordered_map<V, unsigned> DepthMap;
	DepthMap depthMap;
	std::deque<V> q;

	q.push_back(start);

	visited.insert(start);
	std::pair<typename DepthMap::iterator, bool> inserted =
		depthMap.insert(std::make_pair(start, 0));
	assert(inserted.second);

	while (!q.empty()) {
		V u = q.front();
		q.pop_front();
		visited.insert(u);
		typename DepthMap::const_iterator it = depthMap.find(u);
		assert(it != depthMap.end());
		unsigned uDepth = it->second;
		if (uDepth == depth)
			return true;
		if (dir == FORWARD) {
			for (boost::tie(oei, oei_end) = out_edges(u, g);
				oei != oei_end; ++oei) {
				V v = target(*oei, g);
				if (visited.find(v) == visited.end()) {
					visited.insert(v);
					std::pair<typename DepthMap::iterator, bool> inserted =
						depthMap.insert(std::make_pair(v, uDepth+1));
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
						depthMap.insert(std::make_pair(v, uDepth+1));
					assert(inserted.second);
					q.push_back(v);
				}
			}
		}
	}

	return false;
}

/**
 * If the given path has only one possible next/prev vertex in the graph,
 * append/prepend that vertex to the path.
 *
 * @param path the path to extend (a list of vertices)
 * @param dir direction of extension (FORWARD or REVERSE)
 * @param g the graph to use for traversal
 * @param trimLen ignore neighbour vertices with branches
 * shorter than this length [0]
 * @return PathExtensionResult: NO_EXTENSION, HIT_BRANCHING_POINT, or EXTENDED
 */
template <class BidirectionalGraph>
static inline SingleExtensionResult extendPathBySingleVertex(
	Path<typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor>& path,
	Direction dir, const BidirectionalGraph& g, unsigned trimLen = 0)
{
	typedef BidirectionalGraph G;
	typedef boost::graph_traits<G> graph_traits;
	typedef typename graph_traits::vertex_descriptor V;

	typename graph_traits::out_edge_iterator oei, oei_end;
	typename graph_traits::in_edge_iterator iei, iei_end;

	V& u = (dir == FORWARD) ? path.back() : path.front();

	if (dir == FORWARD) {
		if (out_degree(u, g) == 0) {
			return SE_DEAD_END;
		} else if (out_degree(u, g) == 1) {
			const V& v = target(*(out_edges(u, g).first), g);
			path.push_back(v);
			return SE_EXTENDED;
		} else {
			V vNext;
			unsigned branchCount = 0;
			for (boost::tie(oei, oei_end) = out_edges(u, g);
				oei != oei_end; ++oei) {
				const V& v = target(*oei, g);
				if (lookAhead(v, dir, trimLen, g)) {
					branchCount++;
					vNext = v;
				}
			}
			if (branchCount == 0) {
				return SE_DEAD_END;
			} else if (branchCount == 1) {
				path.push_back(vNext);
				return SE_EXTENDED;
			} else {
				assert(branchCount > 1);
				return SE_BRANCHING_POINT;
			}
		}
	} else {
		assert(dir == REVERSE);
		if (in_degree(u, g) == 0) {
			return SE_DEAD_END;
		} else if (in_degree(u, g) == 1) {
			const V& v = source(*(in_edges(u, g).first), g);
			path.push_front(v);
			return SE_EXTENDED;
		} else {
			V vNext;
			unsigned branchCount = 0;
			for (boost::tie(iei, iei_end) = in_edges(u, g);
				iei != iei_end; ++iei) {
				const V& v = source(*iei, g);
				if (lookAhead(v, dir, trimLen, g)) {
					branchCount++;
					vNext = v;
				}
			}
			if (branchCount == 0) {
				return SE_DEAD_END;
			} else if (branchCount == 1) {
				path.push_front(vNext);
				return SE_EXTENDED;
			} else {
				assert(branchCount > 1);
				return SE_BRANCHING_POINT;
			}
		}
	}
}

/**
 * Extend a path up to the next branching point in the graph.
 *
 * @param path path to extend (modified by this function)
 * @param dir direction to extend path (FORWARD or REVERSE)
 * @param g graph in which to perform the extension
 * @param trimLen ignore branches less than this length when
 * detecting branch points [0]
 * @return PathExtensionResult: NO_EXTENSION, HIT_BRANCHING_POINT,
 * or EXTENDED.
 */
template <class BidirectionalGraph>
PathExtensionResult extendPath(
	Path<typename boost::graph_traits<BidirectionalGraph>::vertex_descriptor>& path,
	Direction dir, const BidirectionalGraph& g, unsigned trimLen = 0)
{
	typedef BidirectionalGraph G;
	typedef boost::graph_traits<G> graph_traits;
	typedef typename graph_traits::vertex_descriptor V;
	typename graph_traits::out_edge_iterator oei, oei_end;
	typename graph_traits::in_edge_iterator iei, iei_end;

	assert(path.size() > 0);
	size_t origPathLen = path.size();

	/* track visited nodes to avoid infinite traversal of cycles */
	unordered_set<V> visited;
	visited.insert(path.begin(), path.end());

	SingleExtensionResult result;
	bool detectedCycle = false;
	do {
		result = extendPathBySingleVertex(path, dir, g, trimLen);
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
	} while (result == SE_EXTENDED && !detectedCycle);

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
		} else {
			assert(result == SE_BRANCHING_POINT);
			return EXTENDED_TO_BRANCHING_POINT;
		}
	} else {
		assert(path.size() == origPathLen);
		if (detectedCycle) {
			return CYCLE;
		} else if (result == SE_DEAD_END) {
			return DEAD_END;
		} else {
			assert(result == SE_BRANCHING_POINT);
			return BRANCHING_POINT;
		}
	}
}

#endif
