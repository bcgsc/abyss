#ifndef CONTIGPROPERTIES_H
#define CONTIGPROPERTIES_H 1

#include "ContigNode.h"
#include "Graph/Options.h"
#include "Graph/Properties.h"
#include "IOUtil.h"
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <cassert>
#include <iostream>

using boost::graph_traits;

/** Contig length property. */
struct Length
{
	unsigned length;

	bool operator==(const Length& o) const
	{
		return length == o.length;
	}

	Length& operator+=(const Length& o)
	{
		length += o.length;
		return *this;
	}

	template <typename T>
	Length& operator+=(const T& o)
	{
		assert((int)length + (int)o.distance >= 0);
		length += o.distance;
		return *this;
	}

	friend std::ostream& operator<<(std::ostream& out,
			const Length& o)
	{
		return out << "l=" << o.length;
	}

	friend std::istream& operator>>(std::istream& in, Length& o)
	{
		if (in >> std::ws && in.peek() == 'l')
			return in >> expect("l =") >> o.length;
		else
			return in >> o.length;
	}
};

static inline
void put(vertex_length_t, Length& vp, unsigned length)
{
	vp.length = length;
}

static inline
void put(vertex_coverage_t, Length&, unsigned)
{
}

/** The length and coverage of a contig. */
struct ContigProperties {
	unsigned length;
	unsigned coverage;

	ContigProperties() : length(0), coverage(0) { }
	ContigProperties(unsigned length, unsigned coverage)
		: length(length), coverage(coverage) { }

	bool operator==(const ContigProperties& o) const
	{
		return length == o.length && coverage == o.coverage;
	}

	ContigProperties& operator +=(const ContigProperties& o)
	{
		length += o.length;
		coverage += o.coverage;
		return *this;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const ContigProperties& o)
	{
		float coverage = opt::k <= 0 ? 0
			: (float)o.coverage / (o.length - opt::k + 1);
		switch (opt::format) {
		  case ADJ:
			return out << ' ' << o.length << ' ' << o.coverage;
		  case DOT:
			return out << "l=" << o.length << " C=" << o.coverage;
		  case DOT_MEANCOV:
			assert(opt::k > 0);
			return out << "l=" << o.length << " c=" << coverage;
		  case SAM:
			out << "\tLN:" << o.length;
			return opt::k > 0
				? (out << "\tXc:" << coverage)
				: (out << "\tXC:" << o.coverage);
		}
		return out;
	}

	friend std::istream& operator >>(std::istream& in,
			ContigProperties& o)
	{
		if (in >> std::ws && in.peek() == 'l') {
			in >> expect("l =") >> o.length;
			if (in >> std::ws && in.peek() == 'C')
				in >> expect("C =") >> o.coverage;
			return in;
		} else if (in.peek() == 'C') {
			return in >> expect("C=") >> o.coverage
				>> expect(", l=") >> o.length;
		} else
			return in >> o.length >> o.coverage;
	}
};

static inline
void put(vertex_length_t, ContigProperties& vp, unsigned length)
{
	vp.length = length;
}

static inline
void put(vertex_coverage_t, ContigProperties& vp, unsigned coverage)
{
	vp.coverage = coverage;
}

/** The distance between two contigs. */
struct Distance {
	int distance;
	Distance() : distance(-opt::k + 1) { }
	Distance(int d) : distance(d) { }

	bool operator==(const Distance& o) const
	{
		return distance == o.distance;
	}

	bool operator!=(const Distance& o) const
	{
		return distance != o.distance;
	}

	friend std::ostream& operator<<(std::ostream& out,
			const Distance& o)
	{
		return out << "d=" << o.distance;
	}

	friend std::istream& operator>>(std::istream& in, Distance& o)
	{
		in >> expect(" d = ");
		if (in.peek() == '"')
			return in >> expect("\"") >> o.distance >> expect("\"");
		else
			return in >> o.distance;
	}
};

/** Add the specified distance (overlap) to the specified contig
 * length.
 */
template <typename T>
ContigProperties& operator+=(ContigProperties& a, const T& b)
{
	assert((int)a.length + (int)b.distance > 0);
	a.length += b.distance;
	return a;
}

/** Return the distance between two vertices. */
template <typename Graph>
int get(edge_distance_t, const Graph& g,
		typename graph_traits<Graph>::edge_descriptor e)
{
	return g[e].distance;
}

/** Return the edge properties of (u,v) unless either u or v is
 * ambiguous, in which case return a default-constructed instance of
 * the edge properties.
 */
template <typename Graph>
typename edge_bundle_type<Graph>::type
get(edge_bundle_t, const Graph& g, ContigNode u, ContigNode v)
{
	typedef typename graph_traits<Graph>::edge_descriptor
		edge_descriptor;
	if (u.ambiguous() || v.ambiguous()) {
		return typename edge_bundle_type<Graph>::type();
	} else {
		std::pair<edge_descriptor, bool> e = edge(u, v, g);
		if (!e.second)
			std::cerr << "error: no edge "
				<< get(vertex_name, g, u) << " -> "
				<< get(vertex_name, g, v) << '\n';
		assert(e.second);
		return g[e.first];
	}
}

/** Edge weight property map. */
template <typename Graph>
struct EdgeWeightMap {
	typedef typename Graph::edge_descriptor key_type;
	typedef int value_type;
	typedef value_type reference;
	typedef boost::readable_property_map_tag category;

	EdgeWeightMap(const Graph& g) : m_g(g) { }

	reference operator[](const key_type& e) const
	{
		int weight = m_g[e].distance + m_g[target(e, m_g)].length;
		if (weight < 0)
			std::cerr << "error: invalid edge: "
				<< get(edge_name, m_g, e)
				<< " [" << get(edge_bundle, m_g, e) << "]\n";
		assert(weight >= 0);
		return weight;
	}

  private:
	const Graph& m_g;
};

/** Return the edge weight. */
template <typename Graph>
typename EdgeWeightMap<Graph>::value_type
get(edge_weight_t, const Graph& g,
		typename graph_traits<Graph>::edge_descriptor e)
{
	return EdgeWeightMap<Graph>(g)[e];
}

#endif
