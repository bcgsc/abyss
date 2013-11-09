#ifndef ALLPATHS_SEARCH_H_
#define ALLPATHS_SEARCH_H_

#include "Common/Warnings.h"
#include "Common/UnorderedSet.h"
#include "Graph/Path.h"
#include <boost/tuple/tuple.hpp> // for boost::tie
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>

template <class IncidenceGraph>
PathSearchResult allPathsSearch(
	const IncidenceGraph& g,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor start,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor goal,
	std::vector< Path<typename boost::graph_traits<IncidenceGraph>::vertex_descriptor> >& pathsFound)
{
	return allPathsSearch(g, start, goal, NO_LIMIT, NO_LIMIT, NO_LIMIT, pathsFound);
}

template <class IncidenceGraph>
PathSearchResult allPathsSearch(
	const IncidenceGraph& g,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor start,
	typename boost::graph_traits<IncidenceGraph>::vertex_descriptor goal,
	unsigned maxPaths,
	unsigned minDepth,
	unsigned maxDepth,
	std::vector< Path<typename boost::graph_traits<IncidenceGraph>::vertex_descriptor> >& pathsFound)
{
    BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<IncidenceGraph>));
    typedef typename boost::graph_traits<IncidenceGraph>::vertex_descriptor V;
    typedef typename boost::graph_traits<IncidenceGraph>::out_edge_iterator EdgeIter;
	typedef typename std::pair<EdgeIter,EdgeIter> EdgeIterPair;

	unordered_set<V, std::hash<V> > visited;
	Path<V> path;
	std::vector<EdgeIterPair> eiStack;

	path.push_back(start);
	visited.insert(start);
	eiStack.push_back(out_edges(start, g));

	while(!path.empty()) {

		if (path.back() == goal &&
			(minDepth == NO_LIMIT || (path.size() - 1) >= minDepth)) {
			if (maxPaths != NO_LIMIT && pathsFound.size() >= maxPaths)
				return TOO_MANY_PATHS;
			pathsFound.push_back(path);
		}

		// find next unvisited node and append to path
		while(!path.empty()) {
			if ((maxDepth != NO_LIMIT && (path.size() - 1) >= maxDepth) ||
				eiStack.back().first == eiStack.back().second)
			{
				visited.erase(path.back());
				path.pop_back();
				eiStack.pop_back();
				assert(path.empty() == eiStack.empty());
				if (!path.empty())
					eiStack.back().first++;

			} else {
				V v = target(*(eiStack.back().first), g);
				if (visited.find(v) != visited.end()) {
					eiStack.back().first++;
				} else {
					path.push_back(v);
					eiStack.push_back(out_edges(v, g));
					visited.insert(v);
					break;
				}
			} // if ei != ei_end
		} // while !path.empty()

	} // while !path.empty()

	if (pathsFound.empty())
		return NO_PATH;
	else
		return FOUND_PATH;
}

#endif
