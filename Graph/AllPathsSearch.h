#ifndef ALLPATHS_SEARCH_H_
#define ALLPATHS_SEARCH_H_

#include "Common/UnorderedSet.h"
#include "Graph/Path.h"
#include <boost/tuple/tuple.hpp> // for boost::tie
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <vector>

/** result of exhaustive path search from vertex A to vertex B */
template <class VertexT>
struct AllPathsSearchResult
{
public:

	/** code indicating type of success/failure for search */
	PathSearchResult resultCode;
	/** number of edges traversed during search */
	unsigned cost;
	/** all paths from vertex A to vertex B (if successful) */
	std::vector< Path<VertexT> > paths;

	/** constructor */
	AllPathsSearchResult() : resultCode(NO_PATH), cost(0) {}
};

template <class IncidenceGraph>
AllPathsSearchResult<typename boost::graph_traits<IncidenceGraph>::vertex_descriptor>
allPathsSearch(
	const IncidenceGraph& g,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor start,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor goal,
	unsigned maxPaths, unsigned minDepth, unsigned maxDepth, unsigned maxCost)
{
    typedef typename boost::graph_traits<IncidenceGraph>::vertex_descriptor V;
    typedef typename boost::graph_traits<IncidenceGraph>::out_edge_iterator EdgeIter;
	typedef typename std::pair<EdgeIter,EdgeIter> EdgeIterPair;

	unordered_set<V, hash<V> > visited;
	Path<V> path;
	std::vector<EdgeIterPair> eiStack;
	unordered_set<V> cycleVertices;

	path.push_back(start);
	visited.insert(start);
	eiStack.push_back(out_edges(start, g));

	AllPathsSearchResult<V> result;

	while(!path.empty() && result.cost <= maxCost) {

		if (path.back() == goal &&
			(minDepth == NO_LIMIT || (path.size() - 1) >= minDepth)) {
			if (maxPaths != NO_LIMIT && result.paths.size() >= maxPaths) {
				result.resultCode = TOO_MANY_PATHS;
				return result;
			}
			if (!cycleVertices.empty()) {
				result.resultCode = PATH_CONTAINS_CYCLE;
				return result;
			}
			result.paths.push_back(path);
		}

		// find next unvisited node and append to path
		while(!path.empty()) {
			if ((maxDepth != NO_LIMIT && (path.size() - 1) >= maxDepth) ||
				eiStack.back().first == eiStack.back().second)
			{
				visited.erase(path.back());
				cycleVertices.erase(path.back());
				path.pop_back();
				eiStack.pop_back();
				assert(path.empty() == eiStack.empty());
				if (!path.empty())
					eiStack.back().first++;

			} else {
				V v = target(*(eiStack.back().first), g);
				if (visited.find(v) != visited.end()) {
					cycleVertices.insert(v);
					eiStack.back().first++;
				} else {
					path.push_back(v);
					eiStack.push_back(out_edges(v, g));
					visited.insert(v);
					result.cost++;
					break;
				}
			} // if ei != ei_end
		} // while !path.empty()

	} // while !path.empty()

	if (result.cost > maxCost)
		result.resultCode = MAX_COST_EXCEEDED;
	else if (result.paths.empty())
		result.resultCode = NO_PATH;
	else
		result.resultCode = FOUND_PATH;

	return result;
}

template <class IncidenceGraph>
AllPathsSearchResult<typename boost::graph_traits<IncidenceGraph>::vertex_descriptor>
allPathsSearch(
	const IncidenceGraph& g,
		typename boost::graph_traits<IncidenceGraph>::vertex_descriptor start,
		typename boost::graph_traits<IncidenceGraph>::vertex_descriptor goal)
{
	return allPathsSearch(g, start, goal, NO_LIMIT, NO_LIMIT, NO_LIMIT, NO_LIMIT);
}

#endif
