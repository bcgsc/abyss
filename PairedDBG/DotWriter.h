#ifndef PAIREDDBG_DOTWRITER_H
#define PAIREDDBG_DOTWRITER_H 1

/** Written by Shaun Jackman <sjackman@bcgsc.ca>. */

#include "Graph/ContigGraphAlgorithms.h"
#include "SequenceCollection.h"
#include <cassert>
#include <ostream>

class DotWriter
{
public:
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
	typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
	typedef graph_traits<Graph>::adjacency_iterator adjacency_iterator;

private:

/** Write out the specified contig. */
static
void writeContig(std::ostream& out, const Graph& g, const vertex_descriptor& u)
{
	if (contiguous_in(g, u))
		return;
	unsigned n = 1;
	vertex_descriptor v = u;
	while (contiguous_out(g, v)) {
		n++;
		v = *adjacent_vertices(v, g).first;
	}
	out << '"' << get(vertex_name, g, u) << "\" -> \""
		<< get(vertex_name, g, v) << '"';
	if (n > 2)
		out << " [label=" << n << ']';
	out << '\n';
}

/** Write out the contigs that split at the specified sequence. */
static
void writeEdges(std::ostream& out, const Graph& g, const vertex_descriptor& u)
{
	unsigned outdeg = out_degree(u, g);
	if (outdeg == 0)
		return;
	out << '"' << get(vertex_name, g, u) << "\" ->";
	if (outdeg > 1)
		out << " {";
	std::pair<adjacency_iterator, adjacency_iterator>
		adj = adjacent_vertices(u, g);
	for (adjacency_iterator v = adj.first; v != adj.second; ++v)
		out << " \"" << get(vertex_name, g, *v) << '"';
	if (outdeg > 1)
		out << " }";
	out << '\n';
}

/** Write out a dot graph around the specified sequence. */
static
void write_vertex(std::ostream& out, const Graph& g, const vertex_descriptor& u)
{
	if (contiguous_out(g, u))
		writeContig(out, g, u);
	else
		writeEdges(out, g, u);
}

public:

/** Write out a dot graph for the specified collection. */
static
void write(std::ostream& out, const Graph& g)
{
	out << "digraph g {\n";
	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u) {
		if (get(vertex_removed, g, *u))
			continue;
		write_vertex(out, g, *u);
	}
	out << "}" << std::endl;
}

};

#endif
