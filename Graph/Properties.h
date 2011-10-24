#ifndef GRAPH_PROPERTIES_H
#define GRAPH_PROPERTIES_H 1

#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

/** The coverage of a vertex. */
enum vertex_coverage_t { vertex_coverage };

/** The length of a vertex. */
enum vertex_length_t { vertex_length };

/** A property indicating that this vertex has been removed. */
enum vertex_removed_t { vertex_removed };

/** The orientation of a vertex. */
enum vertex_sense_t { vertex_sense };

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
	BOOST_INSTALL_PROPERTY(vertex, coverage);
	BOOST_INSTALL_PROPERTY(vertex, length);
	BOOST_INSTALL_PROPERTY(vertex, removed);
	BOOST_INSTALL_PROPERTY(vertex, sense);
}

#endif
