#ifndef CONSTRAINEDSEARCH_H
#define CONSTRAINEDSEARCH_H 1

#include "ContigGraph.h"
#include "ContigPath.h"
#include "ContigProperties.h"
#include "DirectedGraph.h"
#include <cassert>
#include <istream>
#include <utility>
#include <vector>

namespace opt {
	extern unsigned maxCost;

	/** Abort the search after visiting maxPaths solutions. */
	static const unsigned maxPaths = 200;
}

typedef ContigGraph<DirectedGraph<ContigProperties, Distance> > Graph;
typedef std::pair<ContigNode, int> Constraint;
typedef std::vector<Constraint> Constraints;
typedef std::vector<ContigPath> ContigPaths;

bool constrainedSearch(const Graph& g,
		graph_traits<Graph>::vertex_descriptor origin,
		Constraints& constraints,
		ContigPaths& superPaths, unsigned& compCost);

#endif
