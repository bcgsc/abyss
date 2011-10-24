#ifndef DEPTHFIRSTSEARCH_H
#define DEPTHFIRSTSEARCH_H 1

#include <boost/graph/depth_first_search.hpp>

using boost::graph_traits;

/**
 * Perform a depth-first search starting first with vertices with
 * deg-(u) = 0 and then visiting any remaining vertices.
 */
template <class Graph, class Visitor, class ColorMap>
void depthFirstSearch(const Graph& g, Visitor vis, ColorMap color)
{
	using boost::color_traits;
	using boost::property_traits;
	using boost::tie;

	typedef graph_traits<Graph> GTraits;
	typedef typename GTraits::vertex_descriptor V;
	typedef typename GTraits::vertex_iterator Vit;
	typedef typename property_traits<ColorMap>::value_type ColorValue;
	const ColorValue white = color_traits<ColorValue>::white();

	// Initialize the vertices.
	Vit uit, ulast;
	for (tie(uit, ulast) = vertices(g); uit != ulast; ++uit) {
		V u = *uit;
		put(color, u, white);
		vis.initialize_vertex(u, g);
	}

	// Visit vertices with deg-(u) = 0.
	for (tie(uit, ulast) = vertices(g); uit != ulast; ++uit) {
		V u = *uit;
		if (get(color, u) == white && in_degree(u, g) == 0) {
			vis.start_vertex(u, g);
			boost::detail::depth_first_visit_impl(g, u, vis, color,
					boost::detail::nontruth2());
		}
	}

	// Visit the remaining vertices.
	for (tie(uit, ulast) = vertices(g); uit != ulast; ++uit) {
		V u = *uit;
		if (get(color, u) == white) {
			vis.start_vertex(u, g);
			boost::detail::depth_first_visit_impl(g, u, vis, color,
					boost::detail::nontruth2());
		}
	}
}

#endif
