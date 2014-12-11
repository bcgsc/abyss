#ifndef CONSTRAINEDSEARCH_H
#define CONSTRAINEDSEARCH_H 1

#include "Common/ContigPath.h"
#include "Common/ContigProperties.h"
#include "Graph/ContigGraph.h"
#include "Graph/DirectedGraph.h"
#include "Graph/Properties.h"
#include <algorithm>
#include <climits> // for INT_MIN
#include <cassert>
#include <istream>
#include <utility>
#include <vector>

namespace opt {
	unsigned maxCost = 100000;

	/** Abort the search after visiting maxPaths solutions. */
	const unsigned maxPaths = 200;
}

typedef std::pair<ContigNode, int> Constraint;
typedef std::vector<Constraint> Constraints;
typedef std::vector<ContigPath> ContigPaths;

/** Compare the distance of two constraints. */
static inline bool compareDistance(
		const Constraint& a, const Constraint& b)
{
	return a.second < b.second;
}

/** Compare the ID of a constraint. */
static inline bool compareID(const Constraint& constraint,
		const ContigNode& key)
{
	return constraint.first < key;
}

/** Find a constraint by ID. */
static inline Constraints::iterator findConstraint(
		Constraints& constraints,
		const ContigNode& key)
{
	Constraints::iterator it = lower_bound(
			constraints.begin(), constraints.end(),
			key, compareID);
	return it != constraints.end()
		&& it->first == key ? it : constraints.end();
}

/** Find paths through the graph that satisfy the constraints.
 * @return false if the search exited early
 */
template <typename Graph, typename vertex_descriptor>
bool constrainedSearch(const Graph& g,
		vertex_descriptor u,
		Constraints& constraints,
		Constraints::const_iterator nextConstraint,
		unsigned satisfied,
		ContigPath& path, ContigPaths& solutions,
		int distance, unsigned& visitedCount)
{
	typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;

	assert(satisfied < constraints.size());
	static const int SATISFIED = INT_MAX;
	if (!path.empty()) {
		vertex_descriptor v = path.back();
		Constraints::iterator it = findConstraint(constraints, v);
		if (it != constraints.end() && it->second != SATISFIED) {
			if (distance > it->second)
				return true; // This constraint cannot be met.

			if (++satisfied == constraints.size()) {
				// All the constraints have been satisfied.
				solutions.push_back(path);
				return solutions.size() <= opt::maxPaths;
			}
			// This constraint has been satisfied.
			int constraint = it->second;
			it->second = SATISFIED;
			if (!constrainedSearch(g, u, constraints,
						nextConstraint, satisfied, path, solutions,
						distance, visitedCount))
				return false;
			it->second = constraint;
			return true;
		}

		if (++visitedCount >= opt::maxCost)
			return false; // Too complex.

		// Check that the next constraint has not been violated.
		while (distance > nextConstraint->second
				&& findConstraint(constraints,
					nextConstraint->first)->second == SATISFIED)
			++nextConstraint; // This constraint is satisfied.
		if (distance > nextConstraint->second)
			return true; // This constraint cannot be met.

		distance += g[v].length;
		u = v;
	}

	path.push_back(vertex_descriptor());
	std::pair<out_edge_iterator, out_edge_iterator> adj = g.out_edges(u);
	for (out_edge_iterator it = adj.first; it != adj.second; ++it) {
		path.back() = target(*it, g);
		if (!constrainedSearch(g, u, constraints,
					nextConstraint, satisfied, path, solutions,
					distance + g[*it].distance, visitedCount))
			return false;
	}
	assert(!path.empty());
	path.pop_back();
	return true;
}

/** Find paths through the graph that satisfy the constraints.
 * @return false if the search exited early
 */
template <typename Graph, typename vertex_descriptor>
bool constrainedSearch(const Graph& g,
		vertex_descriptor v,
		Constraints& constraints, ContigPaths& paths,
		unsigned& cost)
{
    if (constraints.empty())
            return false;

	// Sort the constraints by ID.
	sort(constraints.begin(), constraints.end());

	// Sort the constraints by distance.
	Constraints queue(constraints);
	sort(queue.begin(), queue.end(), compareDistance);

	ContigPath path;
	constrainedSearch(g, v, constraints, queue.begin(), 0,
			path, paths, 0, cost);
	return cost >= opt::maxCost ? false : !paths.empty();
}


#endif
