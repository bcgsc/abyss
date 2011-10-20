#ifndef GRAPH_PROPERTIES_H
#define GRAPH_PROPERTIES_H 1

#include "config.h"
#if HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP
# include <boost/graph/properties.hpp>
# include <boost/property_map/property_map.hpp>
#endif

/** The length of a vertex. */
enum vertex_length_t { vertex_length };

/** A property indicating that this vertex has been removed. */
enum vertex_removed_t { vertex_removed };

/** The orientation of a vertex. */
enum vertex_sense_t { vertex_sense };

#if HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP
using boost::readable_property_map_tag;

using boost::edge_bundle;
using boost::edge_bundle_t;
using boost::edge_weight;
using boost::edge_weight_t;
using boost::no_property;
using boost::vertex_bundle;
using boost::vertex_bundle_t;
using boost::vertex_index;
using boost::vertex_index_t;

namespace boost {
	BOOST_INSTALL_PROPERTY(vertex, length);
	BOOST_INSTALL_PROPERTY(vertex, removed);
	BOOST_INSTALL_PROPERTY(vertex, sense);
}

#else // HAVE_BOOST_GRAPH_GRAPH_TRAITS_HPP

/** Readable property map. */
enum readable_property_map_tag { };

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
