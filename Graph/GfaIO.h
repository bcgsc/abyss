#ifndef GFAIO_H
#define GFAIO_H 1

#include "Common/IOUtil.h"
#include "Graph/Properties.h"
#include <boost/graph/graph_traits.hpp>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

using boost::graph_traits;

struct Distance;
struct DistanceEst;

/** Write a graph in GFA format. */
template <typename Graph>
std::ostream& write_gfa1(std::ostream& out, Graph& g)
{
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename graph_traits<Graph>::edge_iterator Eit;
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::vertex_iterator Vit;
	typedef typename vertex_bundle_type<Graph>::type VP;

	out << "H\tVN:Z:1.0\n";
	assert(out);

	std::pair<Vit, Vit> vrange = vertices(g);
	for (Vit uit = vrange.first; uit != vrange.second; ++uit, ++uit) {
		V u = *uit;
		if (get(vertex_removed, g, u))
			continue;
		const VP& vp = g[u];
		out << "S\t" << get(vertex_contig_name, g, u)
			<< "\t*\tLN:i:" << vp.length;
		if (vp.coverage > 0)
			out << "\tKC:i:" << vp.coverage;
		out << '\n';
	}

	std::pair<Eit, Eit> erange = edges(g);
	for (Eit eit = erange.first; eit != erange.second; ++eit) {
		E e = *eit;
		V u = source(e, g);
		V v = target(e, g);
		if (get(vertex_removed, g, u))
			continue;

		// Output only the canonical edge.
		if (u > get(vertex_complement, g, v))
			continue;

		int distance = g[e].distance;
		out << 'L'
			<< '\t' << get(vertex_contig_name, g, u)
			<< '\t' << (get(vertex_sense, g, u) ? '-' : '+')
			<< '\t' << get(vertex_contig_name, g, v)
			<< '\t' << (get(vertex_sense, g, v) ? '-' : '+');
		if (distance <= 0)
			out << '\t' << -distance << "M\n";
		else
			out << "\t*\n";
		assert(out);
	}
	return out;
}

/** Write a GFA 2 overlap edge. */
template <typename Graph>
std::ostream& write_gfa2_edge(std::ostream& out, const Graph& g,
	typename graph_traits<Graph>::edge_iterator eit)
{
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename graph_traits<Graph>::vertex_descriptor V;

	E e = *eit;
	V u = source(e, g);
	V v = target(e, g);

	int distance = get(edge_bundle, g, eit).distance;
	assert(distance <= 0);
	unsigned overlap = -distance;
	unsigned ulen = g[u].length;
	unsigned vlen = g[v].length;
	bool usense = get(vertex_sense, g, u);
	bool vsense = get(vertex_sense, g, v);
	unsigned ustart = usense ? 0 : ulen - overlap;
	unsigned uend = usense ? overlap : ulen;
	unsigned vstart = !vsense ? 0 : vlen - overlap;
	unsigned vend = !vsense ? overlap : vlen;
	out << "E\t*"
		<< '\t' << get(vertex_name, g, u)
		<< '\t' << get(vertex_name, g, v);
	out << '\t' << ustart;
	if (ustart == ulen)
		out << '$';
	out << '\t' << uend;
	if (uend == ulen)
		out << '$';
	out << '\t' << vstart;
	if (vstart == vlen)
		out << '$';
	out << '\t' << vend;
	if (vend == vlen)
		out << '$';
	out << '\t' << overlap << 'M' << "\n";
	assert(out);
	return out;
}

/** Write GFA 2 gap edge. */
template <typename Graph>
std::ostream& write_gfa2_gap(std::ostream& out, const Graph& g,
	typename graph_traits<Graph>::edge_iterator eit)
{
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename graph_traits<Graph>::vertex_descriptor V;

	E e = *eit;
	V u = source(e, g);
	V v = target(e, g);

	return out << "G\t*"
		<< '\t' << get(vertex_name, g, u)
		<< '\t' << get(vertex_name, g, v)
		<< '\t' << get(edge_bundle, g, eit)
		<< '\n';
}

/** Write GFA 2 overlap edges. */
template <typename Graph>
std::ostream& write_gfa2_edges(std::ostream& out, const Graph& g,
	const Distance*)
{
	typedef typename graph_traits<Graph>::edge_iterator Eit;
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename graph_traits<Graph>::vertex_descriptor V;

	std::pair<Eit, Eit> erange = edges(g);
	for (Eit eit = erange.first; eit != erange.second; ++eit) {
		E e = *eit;
		V u = source(e, g);
		V v = target(e, g);
		if (get(vertex_removed, g, u))
			continue;

		// Output only the canonical edge.
		if (u > get(vertex_complement, g, v))
			continue;
	
		assert(!get(vertex_removed, g, v));
		write_gfa2_edge(out, g, eit);
	}
	return out;
}

/** Write GFA 2 overlap and gap edges. */
template <typename Graph>
std::ostream& write_gfa2_edges(std::ostream& out, const Graph& g,
	const DistanceEst*)
{
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename graph_traits<Graph>::edge_iterator Eit;
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename edge_bundle_type<Graph>::type EP;

	std::pair<Eit, Eit> erange = edges(g);
	for (Eit eit = erange.first; eit != erange.second; ++eit) {
		E e = *eit;
		V u = source(e, g);
		V v = target(e, g);
		if (get(vertex_removed, g, u))
			continue;

		// Output only the canonical edge.
		if (u > get(vertex_complement, g, v))
			continue;

		assert(!get(vertex_removed, g, v));
		const EP& ep = get(edge_bundle, g, eit);
		if (ep.stdDev > 0)
			write_gfa2_gap(out, g, eit);
		else
			write_gfa2_edge(out, g, eit);
	}
	return out;
}

/** Write a graph in GFA 2 format. */
template <typename Graph>
std::ostream& write_gfa2(std::ostream& out, Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename graph_traits<Graph>::vertex_iterator Vit;
	typedef typename vertex_bundle_type<Graph>::type VP;
	typedef typename edge_bundle_type<Graph>::type EP;

	out << "H\tVN:Z:2.0\n";
	assert(out);

	std::pair<Vit, Vit> vrange = vertices(g);
	for (Vit uit = vrange.first; uit != vrange.second; ++uit, ++uit) {
		V u = *uit;
		if (get(vertex_removed, g, u))
			continue;
		const VP& vp = g[u];
		out << "S\t" << get(vertex_contig_name, g, u)
			<< '\t' << vp.length << "\t*";
		if (vp.coverage > 0)
			out << "\tKC:i:" << vp.coverage;
		out << '\n';
	}

	return write_gfa2_edges(out, g, (EP*)NULL);
}

/** Read a graph in GFA format. */
template <typename Graph, typename BetterEP>
std::istream& read_gfa(std::istream& in, Graph& g, BetterEP betterEP)
{
	assert(in);

	typedef typename graph_traits<Graph>::vertex_descriptor V;
	typedef typename vertex_property<Graph>::type VP;
	typedef typename graph_traits<Graph>::edge_descriptor E;
	typedef typename edge_property<Graph>::type EP;

	// Add vertices if this graph is empty.
	bool addVertices = num_vertices(g) == 0;

	while (in && in.peek() != EOF) {
		switch (in.peek()) {
		  case 'H':
			in >> expect("H\tVN:Z:");
			if (in.peek() == '1')
				in >> expect("1.0\n");
			else
				in >> expect("2.0\n");
			assert(in);
			break;

		  case 'S': {
			std::string uname;
			in >> expect("S\t") >> uname >> expect("\t");
			assert(in);

			std::string seq;
			unsigned length = 0;
			if (isdigit(in.peek())) {
				// GFA 2
				in >> length >> seq;
				assert(in);
			} else {
				// GFA 1
				in >> seq;
				assert(in);
				assert(!seq.empty());
				if (seq == "*") {
					in >> expect(" LN:i:") >> length;
					assert(in);
				} else
					length = seq.size();
			}

			unsigned coverage = 0;
			if (in.peek() == '\t' && in.get() == '\t' && in.peek() == 'K') {
				in >> expect("KC:i:") >> coverage;
				assert(in);
			}

			in >> Ignore('\n');
			assert(in);

			if (addVertices) {
				VP vp;
				put(vertex_length, vp, length);
				put(vertex_coverage, vp, coverage);
				V u = add_vertex(vp, g);
				put(vertex_name, g, u, uname);
			} else {
				V u = find_vertex(uname, false, g);
				assert(get(vertex_index, g, u) < num_vertices(g));
				(void)u;
			}
			break;
		  }

		  case 'L': {
			std::string uname, vname;
			char usense, vsense;
			int overlap;
			in >> expect("L\t")
				>> uname >> usense
				>> vname >> vsense >> std::ws;
			if (in.peek() == '*') {
				in.get();
				overlap = -1;
			} else {
				in >> overlap >> expect("M");
			}
			in >> Ignore('\n');
			assert(in);
			assert(!uname.empty());
			assert(!vname.empty());
			assert(usense == '+' || usense == '-');
			assert(vsense == '+' || vsense == '-');

			V u = find_vertex(uname, usense == '-', g);
			V v = find_vertex(vname, vsense == '-', g);
			if (overlap >= 0) {
				int d = -overlap;
				EP ep(d);
				add_edge(u, v, ep, g);
			} else
				add_edge(u, v, g);
			break;
		  }

		  case 'E': {
			std::string ename, uname, vname;
			in >> expect("E\t") >> ename >> uname >> vname;
			assert(in);
			unsigned ustart, uend, vstart, vend;
			in >> ustart >> Skip('$')
				>> uend >> Skip('$')
				>> vstart >> Skip('$')
				>> vend >> Skip('$')
				>> Ignore('\n');
			assert(in);
			unsigned ulength = uend - ustart;
			unsigned vlength = vend - vstart;
			if (ulength != vlength) {
				std::cerr << "error: alignment contains gaps: " << ename << '\t' << uname << '\t' << vname << '\n';
				exit(EXIT_FAILURE);
			}
			V u = find_vertex(uname, g);
			V v = find_vertex(vname, g);
			int d = -ulength;
			EP ep(d);
			add_edge(u, v, ep, g);
			break;
		  }

		  case 'G': {
			std::string ename, uname, vname;
			in >> expect("G\t") >> ename >> uname >> vname;
			assert(in);
			EP ep;
			in >> ep >> Ignore('\n');
			assert(in);
			V u = find_vertex(uname, g);
			V v = find_vertex(vname, g);
			E e;
			bool found;
			boost::tie(e, found) = edge(u, v, g);
			if (found) {
				// Parallel edge
				EP& ref = g[e];
				ref = betterEP(ref, ep);
			} else
				add_edge(u, v, ep, g);
			break;
		  }

		  case '#': // comment
		  case 'C': // GFA1 containment
		  case 'F': // GFA2 fragment
		  case 'O': // GFA2 ordered path
		  case 'P': // GFA1 path
		  case 'U': // GFA2 unordered set
			in >> Ignore('\n');
			break;

		  default: {
			std::string s;
			in >> s >> Ignore('\n');
			std::cerr << "warning: unknown record type: `" << s << "'\n";
		  }
		}
	}
	assert(in.eof());
	return in;
}

#endif
