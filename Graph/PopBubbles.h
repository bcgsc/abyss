#ifndef POPBUBBLES_H
#define POPBUBBLES_H 1

#include "Common/Options.h"
#include "DepthFirstSearch.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <set>
#include <utility>
#include <vector>

using boost::graph_traits;

/** Record a topological order of the vertices. */
template <typename OutIt>
struct TopoVisitor : public boost::default_dfs_visitor
{
	TopoVisitor(OutIt it) : m_it(it) { }

	template <typename Vertex, typename Graph>
	void finish_vertex(const Vertex& u, Graph&) { *m_it++ = u; }

  private:
	OutIt m_it;
};

/** Record a topological order of the vertices. */
template <typename Graph, typename It>
void topologicalSort(const Graph& g, It it)
{
	using boost::default_color_type;
	using boost::property_map;
	using boost::vector_property_map;
	typedef typename property_map<Graph, vertex_index_t>::type
		VertexIndexMap;
	typedef vector_property_map<default_color_type, VertexIndexMap>
		ColorMap;
	depthFirstSearch(g, TopoVisitor<It>(it),
			ColorMap(num_vertices(g)));
}

/** Return true if the specified sequence of vertices is a bubble. */
template <typename Graph, typename It>
bool isBubble(const Graph& g, It first, It last)
{
	typedef typename graph_traits<Graph>::adjacency_iterator Ait;
	typedef typename graph_traits<Graph>::in_edge_iterator Iit;
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	assert(last - first > 1);
	if (last - first == 2)
		return false; // unambiguous edge
	if (*first == get(vertex_complement, g, last[-1]))
		return false; // palindrome
	std::set<V> targets(first, first + 1);
	for (It it = first; it != last - 1; ++it) {
		std::pair<Ait, Ait> adj = adjacent_vertices(*it, g);
		targets.insert(adj.first, adj.second);
	}
	std::set<V> sources(last - 1, last);
	for (It it = first + 1; it != last; ++it) {
		std::pair<Iit, Iit> adj = in_edges(*it, g);
		for (Iit eit = adj.first; eit != adj.second; ++eit)
			sources.insert(source(*eit, g));
	}
	std::set<V> bubble(first, last);
	return sources == bubble && targets == bubble;
}

typedef std::vector<ContigNode> Bubble;
typedef std::vector<Bubble> Bubbles;

/** Discover bubbles. */
template <typename Graph>
Bubbles discoverBubbles(const Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;

	std::vector<V> topo(num_vertices(g));
	topologicalSort(g, topo.rbegin());

	Bubbles bubbles;
	typedef typename std::vector<V>::const_iterator It;
	for (It first = topo.begin(); first != topo.end(); ++first) {
		int sum = out_degree(*first, g);
		if (sum < 2)
			continue;
		if (opt::verbose > 3)
			std::cerr << "* " << get(vertex_name, g, *first) << '\n';
		for (It it = first + 1; it != topo.end(); ++it) {
			unsigned indeg = in_degree(*it, g);
			unsigned outdeg = out_degree(*it, g);
			sum -= indeg;

			if (opt::verbose > 3)
				std::cerr << get(vertex_name, g, *it)
					<< '\t' << indeg << '\t' << outdeg
					<< '\t' << sum
					<< '\t' << sum + (int)outdeg << '\n';

			if (indeg == 0 || sum < 0)
				break;
			if (sum == 0) {
				It last = it + 1;
				if (isBubble(g, first, last)) {
					if (opt::verbose > 3)
						std::cerr << "good\n";
					bubbles.push_back(Bubble(first, last));
					first = it - 1;
				}
				break;
			}

			if (outdeg == 0)
				break;
			sum += outdeg;
		}
	}
	return bubbles;
}

/** Return the length of the longest path through the bubble. */
template <typename Graph>
int longestPath(const Graph& g, const Bubble& topo)
{
	using boost::tie;
	typedef graph_traits<Graph> GTraits;
	typedef typename GTraits::edge_descriptor E;
	typedef typename GTraits::out_edge_iterator Eit;
	typedef typename GTraits::vertex_descriptor V;

	EdgeWeightMap<Graph> weight(g);
	/* GCC 4.4.6 has a bug that prevents ContigNode being used as the
	 * key of a std::map. Possibly related to
	 * http://gcc.gnu.org/bugzilla/show_bug.cgi?id=39390
	 */
	boost::unordered_map<V, int> distance;
	distance[topo.front()] = 0;
	for (Bubble::const_iterator it = topo.begin();
			it != topo.end(); ++it) {
		V u = *it;
		Eit eit, elast;
		for (tie(eit, elast) = out_edges(u, g); eit != elast; ++eit) {
			E e = *eit;
			V v = target(e, g);
			distance[v] = std::max(
					distance[v], distance[u] + weight[e]);
		}
	}
	V v = topo.back();
	return distance[v] - g[v].length;
}

/** Scaffold over the bubble between vertices u and w.
 * Add an edge (u,w) with the distance property set to the length of
 * the largest branch of the bubble.
 */
template <typename Graph>
void scaffoldBubble(Graph& g, const Bubble& bubble)
{
	typedef graph_traits<Graph> GTraits;
	typedef typename GTraits::vertex_descriptor V;
	assert(bubble.size() > 2);

	V u = bubble.front(), w = bubble.back();
	if (edge(u, w, g).second) {
		// Already scaffolded.
		return;
	}
	assert(isBubble(g, bubble.begin(), bubble.end()));

	add_edge(u, w, std::max(longestPath(g, bubble), 1), g);
}

/** Replace each bubble in the graph with a single edge.
 * Remove the vertices in the bubbles from the graph.
 * @return the vertices that were removed from the graph
 */
template <typename Graph>
std::vector<typename graph_traits<Graph>::vertex_descriptor>
popBubbles(Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef std::vector<V> Vertices;
	Vertices popped;
	Bubbles bubbles = discoverBubbles(g);
	for (Bubbles::const_iterator it = bubbles.begin();
			it != bubbles.end(); ++it) {
		scaffoldBubble(g, *it);
		popped.insert(popped.end(), it->begin() + 1, it->end() - 1);
	}
	for (typename Vertices::const_iterator it = popped.begin();
			it != popped.end(); ++it) {
		V u = *it;
		clear_vertex(u, g);
		remove_vertex(u, g);
	}
	return popped;
}

#endif
