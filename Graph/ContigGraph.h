#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "DirectedGraph.h"
#include <istream>
#include <limits> // for numeric_limits
#include <ostream>
#include <sstream>

/** A contig graph is a directed graph with the property that
 * the edge (u,v) implies the existence of the edge (~v,~u).
 */
template <typename VertexProp = no_property>
class ContigGraph : public DirectedGraph<VertexProp> {
  public:
	typedef DirectedGraph<VertexProp> DG;

	// Vertex types.
	typedef typename DG::vertices_size_type vertices_size_type;
	typedef typename DG::vertex_descriptor vertex_descriptor;
	typedef typename DG::vertex_iterator vertex_iterator;
	typedef typename DG::Vertex Vertex;

	// Edge types.
	typedef typename DG::degree_size_type degree_size_type;
	typedef typename DG::out_edge_iterator out_edge_iterator;

  public:
	/** Construct an empty contig graph. */
	ContigGraph() { }

	/** Construct a contig graph with n vertices. The underlying
	 * directed graph has two vertices for each contig. */
	ContigGraph(vertices_size_type n) : DG(2 * n) { }

	/** Return the in degree of vertex v. */
	degree_size_type in_degree(vertex_descriptor v) const
	{
		return (*this)[~v].out_degree();
	}

	/** Return the in degree of vertex v. */
	degree_size_type in_degree(const Vertex& v) const
	{
		return in_degree(vertex(v));
	}

	/** Remove all out edges from vertex v. */
	void clear_out_edges(vertex_descriptor v)
	{
		const Vertex& vertex = (*this)[v];
		for (out_edge_iterator it = vertex.begin();
				it != vertex.end(); ++it)
			remove_edge(~target(*it), ~v);
		DG::clear_out_edges(v);
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

  private:
	ContigGraph(const ContigGraph&);
};

/** Output a contig adjacency graph. */
template <typename VertexProp>
std::ostream& operator<<(std::ostream& out,
		const ContigGraph<VertexProp>& g)
{
	typedef ContigGraph<VertexProp> G;
	for (typename G::vertex_iterator v = g.begin();
			v != g.end(); ++v) {
		const ContigNode& id = g.vertex(*v);
		if (!id.sense())
			out << g.vertex(*v).id();
	   	out << "\t;";
		for (typename G::out_edge_iterator e = v->begin();
				e != v->end(); ++e)
			out << ' ' << (g.target(*e) ^ id.sense());
		if (id.sense())
			out << '\n';
	}
	return out;
}

/** Read the edges of a vertex. */
template <typename VertexProp>
void readEdges(std::istream& in, LinearNumKey id,
		ContigGraph<VertexProp>& graph)
{
	for (int sense = false; sense <= true; ++sense) {
		std::string s;
		getline(in, s, !sense ? ';' : '\n');
		assert(in.good());
		std::istringstream ss(s);
		for (ContigNode edge; ss >> edge;)
			graph.add_edge(ContigNode(id, sense), edge ^ sense);
		assert(ss.eof());
	}
}

/** Read a contig adjacency graph. */
template <typename VertexProp>
std::istream& operator>>(std::istream& in, ContigGraph<VertexProp>& o)
{
	assert(!g_contigIDs.empty());

	// Create the vertices.
	o.clear();
	ContigGraph<VertexProp>(g_contigIDs.size()).swap(o);

	// Load the edges.
	assert(in);
	for (std::string id; in >> id;) {
		in.ignore(std::numeric_limits<std::streamsize>::max(), ';');
		readEdges(in, stringToID(id), o);
	}
	assert(in.eof());
	return in;
}

void readContigGraph(ContigGraph<>& graph, const std::string& path);

#endif
