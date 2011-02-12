/** Written by Shaun Jackman <sjackman@bcgsc.ca>. */

#include "DotWriter.h"
#include "ContigGraphAlgorithms.h"
#include "SequenceCollection.h"
#include <cassert>
#include <ostream>

using namespace std;

typedef SequenceCollectionHash Graph;
typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
typedef graph_traits<Graph>::adjacency_iterator adjacency_iterator;

/** Write out the specified contig. */
static void write_contig(ostream& out, const Graph& g, const Kmer& u)
{
	if (contiguous_in(g, u))
		return;
	unsigned n = 1;
	Kmer v = u;
	while (contiguous_out(g, v)) {
		n++;
		v = *adjacent_vertices(v, g).first;
	}
	if (n > 1)
		out << u << "->" << v << "[label=\"" << n << "\"]\n";
}

/** Write out the contigs that split at the specified sequence. */
static void write_split(ostream& out, const Graph& g, const Kmer& u)
{
	if (out_degree(u, g) <= 1)
		return;
	out << u << "->{";
	std::pair<adjacency_iterator, adjacency_iterator>
		adj = adjacent_vertices(u, g);
	for (adjacency_iterator v = adj.first; v != adj.second; ++v)
		out << ' ' << *v;
	out << " }\n";
}

/** Write out the contigs that join at the specified sequence. */
static void write_join(ostream& out, const Graph& g, const Kmer& u)
{
	if (in_degree(u, g) <= 1)
		return;
	out << "{";
	std::pair<adjacency_iterator, adjacency_iterator>
		adj = adjacent_vertices(~u, g);
	for (adjacency_iterator v = adj.first; v != adj.second; ++v)
		out << ' ' << ~*v;
	out << " }->" << u << '\n';
}

/** Write out a dot graph around the specified sequence. */
static void write_vertex(ostream& out, const Graph& g, const Kmer& u)
{
	write_split(out, g, u);
	write_join(out, g, u);
	write_contig(out, g, u);
}

/** Write out a dot graph for the specified collection. */
void DotWriter::write(ostream& out, const Graph& g)
{
	out << "digraph g {\n";
	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u) {
		if (get(vertex_removed, g, *u))
			continue;
		write_vertex(out, g, *u);
	}
	out << "}" << endl;
}
