#ifndef GRAPH_PROPERTIES_H
#define GRAPH_PROPERTIES_H 1

#include <istream>
#include <ostream>
#include <boost/graph/properties.hpp>

/** The distance between two vertices. */
enum edge_distance_t { edge_distance };

/** The complementary vertex of a skew-symmetric graph. */
enum vertex_complement_t { vertex_complement };

/** The index of a contig. */
enum vertex_contig_index_t { vertex_contig_index };

/** The name of a contig without an orientation. */
enum vertex_contig_name_t { vertex_contig_name };

/** The coverage of a vertex. */
enum vertex_coverage_t { vertex_coverage };

/** The length of a vertex. */
enum vertex_length_t { vertex_length };

/** A property indicating that this vertex has been removed. */
enum vertex_removed_t { vertex_removed };

/** The orientation of a vertex. */
enum vertex_sense_t { vertex_sense };

using boost::edge_bundle;
using boost::edge_bundle_t;
using boost::edge_name;
using boost::edge_name_t;
using boost::edge_weight;
using boost::edge_weight_t;
using boost::no_property;
using boost::vertex_bundle;
using boost::vertex_bundle_t;
using boost::vertex_index;
using boost::vertex_index_t;
using boost::vertex_name;
using boost::vertex_name_t;

using boost::edge_bundle_type;
using boost::edge_property;
using boost::vertex_bundle_type;
using boost::vertex_property;

namespace boost {
	BOOST_INSTALL_PROPERTY(edge, distance);
	BOOST_INSTALL_PROPERTY(vertex, complement);
	BOOST_INSTALL_PROPERTY(vertex, contig_index);
	BOOST_INSTALL_PROPERTY(vertex, contig_name);
	BOOST_INSTALL_PROPERTY(vertex, coverage);
	BOOST_INSTALL_PROPERTY(vertex, length);
	BOOST_INSTALL_PROPERTY(vertex, removed);
	BOOST_INSTALL_PROPERTY(vertex, sense);
}

/** No property. */
struct NoProperty
{
	NoProperty(...) { }
	bool operator==(const NoProperty&) const { return true; }
	friend std::ostream& operator<<(std::ostream& out,
			const NoProperty&)
	{
		return out;
	}
	friend std::istream& operator>>(std::istream& in, NoProperty&)
	{
		return in;
	}
};

template <typename Tag>
void put(Tag, NoProperty&, unsigned)
{
}

#endif
