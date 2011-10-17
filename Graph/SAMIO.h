#ifndef SAMIO_H
#define SAMIO_H 1

#include "ContigNode.h"
#include "ContigProperties.h" // for edge_distance
#include "Graph.h"
#include <ostream>

/** Output a graph in SAM alignment format.
 * vertex_descriptor must be convertible to a ContigNode
 * vertex_property_type must have members length and coverage
 * edge_property_type must have a member distance
 */
template <typename Graph>
std::ostream& write_sam(std::ostream& out, const Graph& g,
		const std::string& program, const std::string& commandLine)
{
	typedef typename graph_traits<Graph>::vertex_iterator
		vertex_iterator;
	typedef typename vertex_property<Graph>::type
		vertex_property_type;
	typedef typename graph_traits<Graph>::edge_iterator
		edge_iterator;
	typedef typename edge_property<Graph>::type
		edge_property_type;

	out << "@HD\tVN:1.0\n"
		"@PG\tID:" << program << "\tVN:" VERSION "\t"
		"CL:" << commandLine << '\n';

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
		ContigNode u = source(*e, g), v = target(*e, g);
		assert(!get(vertex_removed, g, v));
		int distance = get(edge_distance, g, *e);
		if (get(vertex_removed, g, u) || distance >= 0)
			continue;
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

#endif
