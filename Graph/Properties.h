#ifndef GRAPH_PROPERTIES_H
#define GRAPH_PROPERTIES_H 1

#include "config.h"
#if HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP
# include <boost/graph/properties.hpp>
#endif

/** A property indicating that this vertex has been removed. */
enum vertex_removed_t { vertex_removed };

#if HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP
using boost::edge_bundle;
using boost::edge_bundle_t;
using boost::no_property;
using boost::vertex_bundle;
using boost::vertex_bundle_t;
using boost::vertex_index;
using boost::vertex_index_t;

namespace boost {
	BOOST_INSTALL_PROPERTY(vertex, removed);
}

#else // HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP

/** No properties. */
struct no_property { };

/** A vertex bundle property. */
enum vertex_bundle_t { vertex_bundle };

/** The vertex index property. */
enum vertex_index_t { vertex_index } ;

/** An edge bundle property. */
enum edge_bundle_t { edge_bundle };

#endif // HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP

#endif
