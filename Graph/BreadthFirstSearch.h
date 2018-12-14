#ifndef DEPTHFIRSTSEARCH_H
#define DEPTHFIRSTSEARCH_H 1

#include "Graph/DefaultColorMap.h"
#include <vector>
#include <boost/graph/breadth_first_search.hpp>

using boost::function_requires;
using boost::graph_traits;
using boost::property_traits;
using boost::color_traits;

/** return code from BFS visitor callbacks */
enum BFSVisitorResult { BFS_SUCCESS = 0, BFS_SKIP_ELEMENT, BFS_ABORT_SEARCH };

/**
 * Default BFS visitor that does nothing.  Used as a placeholder when
 * no special visitor behaviour is needed.
 */
template <class G>
class DefaultBFSVisitor
{
public:

	typedef typename boost::graph_traits<G>::vertex_descriptor V;
	typedef typename boost::graph_traits<G>::edge_descriptor E;

	BFSVisitorResult discover_vertex(const V&, const G&) { return BFS_SUCCESS; }
	BFSVisitorResult examine_vertex(const V&, const G&) { return BFS_SUCCESS; }
	BFSVisitorResult finish_vertex(const V&, const G&) { return BFS_SUCCESS; }
	BFSVisitorResult examine_edge(const E&, const G&) { return BFS_SUCCESS; }
	BFSVisitorResult tree_edge(const E&, const G&) { return BFS_SUCCESS; }
	BFSVisitorResult non_tree_edge(const E&, const G&) { return BFS_SUCCESS; }
	BFSVisitorResult gray_target(const E&, const G&) { return BFS_SUCCESS; }
	BFSVisitorResult black_target(const E&, const G&) { return BFS_SUCCESS; }
	void post_processing() {}
};

template <class Graph, class Queue, class ColorMap, class Visitor>
static inline BFSVisitorResult bfsVisitEdge(
	const typename boost::graph_traits<Graph>::edge_descriptor& e,
	bool isOutEdge, const Graph& g, Queue& queue, ColorMap& colorMap,
	Visitor& visitor)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;

	typedef typename property_traits<ColorMap>::value_type ColorValue;
	typedef color_traits<ColorValue> Color;

	BFSVisitorResult result = BFS_SUCCESS;

	/* note: boost edges always go from source to target */
	const V& v = isOutEdge ? target(e, g) : source(e, g);

	result = visitor.examine_edge(e, g);
	if (result != BFS_SUCCESS)
		return result;

	ColorValue color = get(colorMap, v);
	if (color == Color::white()) {

		result = visitor.tree_edge(e, g);
		if (result != BFS_SUCCESS)
			return result;

		result = visitor.discover_vertex(v, g);
		if (result != BFS_SUCCESS)
			return result;

		put(colorMap, v, Color::gray());
		queue.push(v);

	} else {

		result = visitor.non_tree_edge(e, g);
		if (result != BFS_SUCCESS)
			return result;

		if (color == Color::gray()) {
			result = visitor.gray_target(e, g);
		} else {
			result = visitor.black_target(e, g);
		}

	}

	return result;
}

/**
 * An implementation of BFS that allows multiple start nodes and permits
 * terminating the search early.  Visitors can terminate the search by
 * returning the BFS_ABORT_SEARCH result code from their callback routines.
 */
template <class Graph, class VertexSet, class ColorMap, class Visitor>
static inline void breadthFirstSearchImpl(
	const VertexSet& startVertices, const Graph& g, bool undirected,
	ColorMap& colorMap, Visitor& visitor)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor V;
	typedef typename boost::graph_traits<Graph>::edge_descriptor E;
	typedef typename boost::graph_traits<Graph>::in_edge_iterator IEIt;
	typedef typename boost::graph_traits<Graph>::out_edge_iterator OEIt;

	typedef typename VertexSet::const_iterator VertexListConstIt;

	typedef typename property_traits<ColorMap>::value_type ColorValue;
	typedef color_traits<ColorValue> Color;

	/* breadth-first search queue */
	boost::queue<V> queue;

	BFSVisitorResult result;
	IEIt iei, iei_end;
	OEIt oei, oei_end;

	/* push start vertices onto search queue */

	for (VertexListConstIt it = startVertices.begin(); it != startVertices.end(); ++it) {
		ColorValue color = get(colorMap, *it);
		if (color == Color::white()) {
			result = visitor.discover_vertex(*it, g);
			if (result == BFS_SKIP_ELEMENT)
				continue;
			else if (result == BFS_ABORT_SEARCH)
				return;
			put(colorMap, *it, Color::gray());
		}
		/*
		 * note: vertex may already be black if user is reusing colorMap across
		 * multiple searches
		 */
		if (color != Color::black()) {
			queue.push(*it);
		}
	}

	/* do breadth first search */

	while (!queue.empty()) {

		V u = queue.top();
		queue.pop();
		result = visitor.examine_vertex(u, g);
		if (result == BFS_SKIP_ELEMENT)
			continue;
		else if (result == BFS_ABORT_SEARCH)
			return;

		for (boost::tie(oei, oei_end) = out_edges(u, g); oei != oei_end; ++oei) {
			const E& e = *oei;
			result = bfsVisitEdge(e, true, g, queue, colorMap, visitor);
			if (result == BFS_ABORT_SEARCH)
				return;
		}

		if (undirected) {
			for (boost::tie(iei, iei_end) = in_edges(u, g); iei != iei_end; ++iei) {
				const E& e = *iei;
				result = bfsVisitEdge(e, false, g, queue, colorMap, visitor);
				if (result == BFS_ABORT_SEARCH)
					return;
			}
		}

		put(colorMap, u, Color::black());
		result = visitor.finish_vertex(u, g);
		if (result == BFS_ABORT_SEARCH)
			return;
	}

	visitor.post_processing();
}

template <class Graph>
void breadthFirstSearch(const Graph& g, bool undirected,
	const typename graph_traits<Graph>::vertex_descriptor& root)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	std::vector<V> startVertices(1, root);
	DefaultColorMap<Graph> colorMap;
	DefaultBFSVisitor<Graph> visitor;
	breadthFirstSearchImpl(startVertices, g, undirected, colorMap, visitor);
}

template <class Graph>
void breadthFirstSearch(const Graph& g,
	const typename graph_traits<Graph>::vertex_descriptor& root)
{
	breadthFirstSearch(g, false, root);
}

template <class Graph, class Visitor>
void breadthFirstSearch(
	const typename boost::graph_traits<Graph>::vertex_descriptor& root,
	const Graph& g, bool undirected, Visitor& visitor)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	std::vector<V> startVertices(1, root);
	DefaultColorMap<Graph> colorMap;
	breadthFirstSearchImpl(startVertices, g, undirected, colorMap, visitor);
}

template <class Graph, class Visitor>
void breadthFirstSearch(
	const typename boost::graph_traits<Graph>::vertex_descriptor& root,
	const Graph& g, Visitor& visitor)
{
	breadthFirstSearch(root, g, false, visitor);
}

template <class Graph, class ColorMap, class Visitor>
void breadthFirstSearch(
	const typename boost::graph_traits<Graph>::vertex_descriptor& root,
	const Graph& g, ColorMap& colorMap, Visitor& visitor)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	std::vector<V> startVertices(1, root);
	breadthFirstSearchImpl(startVertices, g, false, colorMap, visitor);
}

template <class Graph, class VertexSet, class Visitor>
void breadthFirstSearchMulti(
	const VertexSet& startVertices,
	const Graph& g, Visitor& visitor)
{
	DefaultColorMap<Graph> colorMap;
	breadthFirstSearchImpl(startVertices, g, false, colorMap, visitor);
}

template <class Graph, class VertexSet, class Visitor>
	void breadthFirstSearchMulti(
	const VertexSet& startVertices,
	const Graph& g, bool undirected, Visitor& visitor)
{
	DefaultColorMap<Graph> colorMap;
	breadthFirstSearchImpl(startVertices, g, undirected, colorMap, visitor);
}

template <class IncidenceGraph, class BFSVisitor, class ColorMap>
	void breadthFirstSearch(const IncidenceGraph& g,
	typename graph_traits<IncidenceGraph>::vertex_descriptor root,
	BFSVisitor& visitor)
{
	typedef typename graph_traits<IncidenceGraph>::vertex_descriptor V;
	DefaultColorMap<IncidenceGraph> colorMap;
	boost::queue<V> q;
	breadthFirstSearch(g, root, q, visitor, colorMap);
}

#endif
