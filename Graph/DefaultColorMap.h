#ifndef DEFAULTCOLORMAP_H
#define DEFAULTCOLORMAP_H

#include "UnorderedMap.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/functional/hash/hash.hpp>

template <typename G>
class DefaultColorMap
{
public:

	typedef typename boost::graph_traits<G>::vertex_descriptor key_type;
	typedef typename boost::graph_traits<G>::vertex_descriptor& reference;
	typedef boost::default_color_type value_type;
	typedef boost::read_write_property_map_tag category;

	unordered_map<key_type, value_type> map;
};

namespace boost {
template <typename G>
struct property_traits< DefaultColorMap<G> > {
	typedef typename DefaultColorMap<G>::key_type key_type;
	typedef typename DefaultColorMap<G>::reference reference;
	typedef typename DefaultColorMap<G>::value_type value_type;
	typedef typename DefaultColorMap<G>::category category;
};
}

template <typename G>
typename DefaultColorMap<G>::value_type
get(DefaultColorMap<G>& colorMap, typename DefaultColorMap<G>::key_type key)
{
	typedef typename DefaultColorMap<G>::key_type key_type;
	typedef typename DefaultColorMap<G>::value_type value_type;

	typename unordered_map<key_type, value_type>::iterator i = colorMap.map.find(key);

	if (i != colorMap.map.end())
		return i->second;

	return boost::white_color;
}

template <typename G>
void
put(DefaultColorMap<G>& colorMap,
	typename DefaultColorMap<G>::key_type key,
	typename DefaultColorMap<G>::value_type value)
{
	colorMap.map[key] = value;
}

#endif
