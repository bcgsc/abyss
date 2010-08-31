#ifndef SAMIO_H
#define SAMIO_H 1

#include "ContigNode.h"
#include "Graph.h"
#include <ostream>

/** Output a graph in SAM alignment format.
 * vertex_descriptor must be convertible to a ContigNode
 * vertex_property_type must have members length and coverage
 * edge_property_type must have a member distance
 */
template <typename Graph>
std::ostream& write_sam(std::ostream& out, const Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_iterator
		vertex_iterator;
	typedef typename vertex_property<Graph>::type
		vertex_property_type;
	typedef typename graph_traits<Graph>::edge_iterator
		edge_iterator;
	typedef typename edge_property<Graph>::type
		edge_property_type;

	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u, ++u) {
 		if (get(vertex_removed, g, *u))
			continue;
		const vertex_property_type& vp = g[*u];
		out << "@SQ\tSN:" << ContigID(*u)
			<< "\tLN:" << vp.length
			<< "\tXC:" << vp.coverage << '\n';
	}

	std::pair<edge_iterator, edge_iterator> eit = edges(g);
	for (edge_iterator e = eit.first; e != eit.second; ++e) {
		int distance = g[*e].distance;
		if (distance >= 0)
			continue;
		ContigNode u = source(*e, g), v = target(*e, g);
		unsigned flag = u.sense() == v.sense() ? 0 : 0x10; //FREVERSE
		unsigned alen = -distance;
		unsigned pos = 1 + (u.sense() ? 0 : g[u].length - alen);
		out << ContigID(v) // QNAME
			<< '\t' << flag // FLAG
			<< '\t' << ContigID(u) // RNAME
			<< '\t' << pos // POS
			<< "\t255\t"; // MAPQ
		// CIGAR
		unsigned clip = g[v].length - alen;
		if (u.sense())
			out << clip << 'H' << alen << "M\t";
		else
			out << alen << 'M' << clip << "H\t";
		// MRNM MPOS ISIZE SEQ QUAL
		out << "*\t0\t0\t*\t*\n";
	}
	return out;
}

template <typename Graph>
struct SAMWriter
{
	const Graph& g;
	SAMWriter(const Graph& g) : g(g) { }
	friend std::ostream& operator<<(std::ostream& out,
			const SAMWriter& o)
	{
		return write_sam<Graph>(out, o.g);
	}
};

template <typename Graph>
SAMWriter<Graph> sam_writer(const Graph& g)
{
	return SAMWriter<Graph>(g);
}

#endif
