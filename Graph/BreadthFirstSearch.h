#ifndef DEPTHFIRSTSEARCH_H
#define DEPTHFIRSTSEARCH_H 1

#include "Graph/DefaultColorMap.h"
#include <vector>
#include <boost/graph/breadth_first_search.hpp>

using boost::function_requires;
using boost::graph_traits;
using boost::property_traits;
using boost::color_traits;

// Note: This is an exact copy of the boost-provided breadth_first_search,
// except that the visitor and color map arguments are passed by reference.

template <class IncidenceGraph, class Buffer, class BFSVisitor,
            class ColorMap>
  void breadthFirstSearch
    (const IncidenceGraph& g,
     typename boost::graph_traits<IncidenceGraph>::vertex_descriptor s,
     Buffer& Q, BFSVisitor& vis, ColorMap& color)
  {
    function_requires< boost::IncidenceGraphConcept<IncidenceGraph> >();
    typedef graph_traits<IncidenceGraph> GTraits;
    typedef typename GTraits::vertex_descriptor Vertex;
    typedef typename GTraits::edge_descriptor Edge;
    function_requires< boost::BFSVisitorConcept<BFSVisitor, IncidenceGraph> >();
    function_requires< boost::ReadWritePropertyMapConcept<ColorMap, Vertex> >();
    typedef typename property_traits<ColorMap>::value_type ColorValue;
    typedef color_traits<ColorValue> Color;
    typename GTraits::out_edge_iterator ei, ei_end;

    put(color, s, Color::gray());             vis.discover_vertex(s, g);
    Q.push(s);
    while (! Q.empty()) {
      Vertex u = Q.top(); Q.pop();            vis.examine_vertex(u, g);
      for (boost::tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
        Vertex v = target(*ei, g);            vis.examine_edge(*ei, g);
        ColorValue v_color = get(color, v);
        if (v_color == Color::white()) {      vis.tree_edge(*ei, g);
          put(color, v, Color::gray());       vis.discover_vertex(v, g);
          Q.push(v);
        } else {                              vis.non_tree_edge(*ei, g);
          if (v_color == Color::gray())       vis.gray_target(*ei, g);
          else                                vis.black_target(*ei, g);
        }
      } // end for
      put(color, u, Color::black());          vis.finish_vertex(u, g);
    } // end while
  } // breadth_first_visit

template <class IncidenceGraph>
void breadthFirstSearch(const IncidenceGraph& g,
		typename graph_traits<IncidenceGraph>::vertex_descriptor root)
{
	typedef typename graph_traits<IncidenceGraph>::vertex_descriptor V;
	DefaultColorMap<IncidenceGraph> colorMap;
	typename boost::default_bfs_visitor visitor;
	boost::queue<V> q;
	breadthFirstSearch(g, root, q, visitor, colorMap);
}

template <class IncidenceGraph, class BFSVisitor, class ColorMap>
void breadthFirstSearch(const IncidenceGraph& g,
		typename graph_traits<IncidenceGraph>::vertex_descriptor root,
		BFSVisitor& visitor, ColorMap& colorMap)
{
	typedef typename graph_traits<IncidenceGraph>::vertex_descriptor V;
	boost::queue<V> q;
	breadthFirstSearch(g, root, q, visitor, colorMap);
}

#endif
