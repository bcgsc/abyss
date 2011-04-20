// Modifications to the original Boost source by Shaun Jackman 2010
//====================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Copyright 2003 Bruce Barr
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//====================================================================
#ifndef DEPTHFIRSTSEARCH_H
#define DEPTHFIRSTSEARCH_H 1

#include <boost/graph/depth_first_search.hpp>

/** A depth-first search that originates at vertices with
 * in_degree == 0.
 */
template <class VertexListGraph, class DFSVisitor, class ColorMap>
void
depthFirstSearch(
		const VertexListGraph& g, DFSVisitor vis, ColorMap color)
{
	using boost::DFSVisitorConcept;
	using boost::color_traits;
	using boost::function_requires;
	using boost::implicit_cast;
	using boost::property_traits;
	using boost::tie;

	typedef typename graph_traits<VertexListGraph>::vertex_descriptor
		Vertex;
	function_requires<DFSVisitorConcept<
		DFSVisitor, VertexListGraph> >();
	typedef typename property_traits<ColorMap>::value_type ColorValue;
	typedef color_traits<ColorValue> Color;

	typename graph_traits<VertexListGraph>::vertex_iterator
		ui, ui_end;
	for (tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
		Vertex u = implicit_cast<Vertex>(*ui);
		put(color, u, Color::white());
		vis.initialize_vertex(u, g);
	}

	for (tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
		Vertex u = implicit_cast<Vertex>(*ui);
		ColorValue u_color = get(color, u);
		if (u_color == Color::white()
				&& in_degree(u, g) == 0) {
			vis.start_vertex(u, g);
			boost::detail::depth_first_visit_impl(g, u, vis, color,
					boost::detail::nontruth2());
		}
	}

	for (tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
		Vertex u = implicit_cast<Vertex>(*ui);
		ColorValue u_color = get(color, u);
		if (u_color == Color::white()) {
			vis.start_vertex(u, g);
			boost::detail::depth_first_visit_impl(g, u, vis, color,
					boost::detail::nontruth2());
		}
	}
}

#endif
