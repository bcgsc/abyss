#ifndef SAMIO_H
#define SAMIO_H 1

#include "Common/ContigNode.h"
#include "Graph/Properties.h"
#include <istream>
#include <sstream>
#include <boost/graph/graph_traits.hpp>
#include <ostream>

using boost::graph_traits;

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

	out << "@HD\tVN:1.0\n"
		"@PG\tID:" << program << "\tVN:" VERSION "\t"
		"CL:" << commandLine << '\n';

	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u, ++u) {
 		if (get(vertex_removed, g, *u))
			continue;
		const vertex_property_type& vp = g[*u];
		out << "@SQ\tSN:" << get(vertex_contig_name, g, *u)
			<< "\tLN:" << vp.length;
		if (vp.coverage > 0)
			out << "\tXC:" << vp.coverage;
		out << '\n';
	}

	std::pair<edge_iterator, edge_iterator> eit = edges(g);
	for (edge_iterator e = eit.first; e != eit.second; ++e) {
		ContigNode u = source(*e, g), v = target(*e, g);
		assert(!get(vertex_removed, g, v));
		int distance = get(edge_distance, g, *e);
		if (get(vertex_removed, g, u) || distance > 0)
			continue;
		unsigned flag = u.sense() == v.sense() ? 0 : 0x10; //FREVERSE
		unsigned alen = -distance;
		unsigned pos = 1 + (u.sense() ? 0 : g[u].length - alen);
		out << get(vertex_contig_name, g, v) // QNAME
			<< '\t' << flag // FLAG
			<< '\t' << get(vertex_contig_name, g, u) // RNAME
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
std::istream& read_sam_header(std::istream& in, Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename vertex_property<Graph>::type VP;
	assert(in);
	for (std::string line; in.peek() == '@' && getline(in, line);) {
		std::istringstream ss(line);
		std::string type;
		ss >> type;
		if (type != "@SQ")
			continue;

		std::string s;
		VP vp;
		ss >> expect(" SN:") >> s >> expect(" LN:") >> vp;
		assert(ss);

		V u = add_vertex(vp, g);
		put(vertex_name, g, u, s);
	}
	if (g.num_vertices() == 0) {
		std::cerr << "error: no @SQ records in the SAM header\n";
		exit(EXIT_FAILURE);
	}
	return in;
}

#endif
