#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "DirectedGraph.h"
#include <cassert>
#include <istream>
#include <limits> // for numeric_limits
#include <ostream>
#include <sstream>
#include <utility>

template <typename VertexProp = no_property> class ContigGraph;

template <typename VertexProp>
std::ostream& operator<<(std::ostream& out,
		const ContigGraph<VertexProp>& g);

/** A contig graph is a directed graph with the property that
 * the edge (u,v) implies the existence of the edge (~v,~u).
 */
template <typename VertexProp>
class ContigGraph : public DirectedGraph<VertexProp> {
  public:
	typedef DirectedGraph<VertexProp> DG;

	// Vertex types.
	typedef typename DG::vertices_size_type vertices_size_type;
	typedef typename DG::vertex_descriptor vertex_descriptor;
	typedef typename DG::vertex_iterator vertex_iterator;

	// Edge types.
	typedef typename DG::degree_size_type degree_size_type;
	typedef typename DG::edge_descriptor edge_descriptor;
	typedef typename DG::adjacency_iterator adjacency_iterator;

  public:
	/** Construct an empty contig graph. */
	ContigGraph() { }

	/** Construct a contig graph with n vertices. The underlying
	 * directed graph has two vertices for each contig. */
	ContigGraph(vertices_size_type n) : DG(2 * n) { }

	/** Return the in degree of vertex v. */
	degree_size_type in_degree(vertex_descriptor v) const
	{
		return DG::out_degree(~v);
	}

	/** Remove all out edges from vertex u. */
	void clear_out_edges(vertex_descriptor u)
	{
		std::pair<adjacency_iterator, adjacency_iterator>
			adj = DG::adjacent_vertices(u);
		for (adjacency_iterator v = adj.first; v != adj.second; ++v)
			DG::remove_edge(~*v, ~u);
		DG::clear_out_edges(u);
	}

	/** Remove all in edges from vertex v. */
	void clear_in_edges(vertex_descriptor v)
	{
		clear_out_edges(~v);
	}

	/** Remove all edges to and from vertex v. */
	void clear_vertex(vertex_descriptor v)
	{
		clear_out_edges(v);
		clear_in_edges(v);
	}

	/** Add a vertex to this graph. */
	vertex_descriptor add_vertex(
				const VertexProp& data = VertexProp())
	{
		vertex_descriptor v = DG::add_vertex(data);
		DG::add_vertex(data);
		return v;
	}

	/** Remove vertex v from this graph. It is assumed that there
	 * are no edges to or from vertex v. It is best to call
	 * clear_vertex before remove_vertex.
	 */
	void remove_vertex(vertex_descriptor v)
	{
		DG::remove_vertex(v);
		DG::remove_vertex(~v);
	}

	/** Add edge (u,v) to this graph. */
	std::pair<edge_descriptor, bool>
	add_edge(vertex_descriptor u, vertex_descriptor v)
	{
		std::pair<edge_descriptor, bool> e = DG::add_edge(u, v);
		DG::add_edge(~v, ~u);
		return e;
	}

	friend std::ostream& operator<< <>(std::ostream& out,
			const ContigGraph& g);

  private:
	ContigGraph(const ContigGraph&);
};

/** Output a contig adjacency graph. */
template <typename VertexProp>
std::ostream& operator<<(std::ostream& out,
		const ContigGraph<VertexProp>& g)
{
	typedef ContigGraph<VertexProp> G;
	typedef typename G::vertex_iterator vertex_iterator;
	typedef typename G::adjacency_iterator adjacency_iterator;
	std::pair<vertex_iterator, vertex_iterator> vit = g.vertices();
	for (vertex_iterator u = vit.first; u != vit.second; ++u) {
		if (g.is_removed(*u))
			continue;
		bool sense = ContigNode(*u).sense();
		if (!sense)
			out << ContigID(*u) << g[*u];
		out << "\t;";
		std::pair<adjacency_iterator, adjacency_iterator>
			adj = g.adjacent_vertices(*u);
		for (adjacency_iterator v = adj.first; v != adj.second; ++v)
			out << ' ' << (*v ^ sense);
		if (sense)
			out << '\n';
	}
	return out;
}

/** Read a contig adjacency graph. */
template <typename VertexProp>
std::istream& operator>>(std::istream& in, ContigGraph<VertexProp>& g)
{
	typedef typename ContigGraph<VertexProp>::DG DG;
	typedef typename
		ContigGraph<VertexProp>::vertex_descriptor vertex_descriptor;

	// Read the vertex properties.
	g.clear();
	assert(in);
	ContigID id(-1);
	VertexProp prop;
	while (in >> id >> prop) {
		in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		vertex_descriptor v = g.add_vertex(prop);
		assert(v == ContigNode(id, false));
	}
	assert(in.eof());
	assert(g.num_vertices() > 0);
	ContigID::lock();

	// Read the edges.
	in.clear();
	in.seekg(std::ios_base::beg);
	assert(in);
	while (in >> id) {
		in.ignore(std::numeric_limits<std::streamsize>::max(), ';');
		for (int sense = false; sense <= true; ++sense) {
			std::string s;
			getline(in, s, !sense ? ';' : '\n');
			assert(in.good());
			std::istringstream ss(s);
			for (ContigNode edge; ss >> edge;)
				g.DG::add_edge(ContigNode(id, sense), edge ^ sense);
			assert(ss.eof());
		}
	}
	assert(in.eof());
	return in;
}

#endif
