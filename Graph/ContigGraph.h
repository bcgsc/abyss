#ifndef CONTIGGRAPH_H
#define CONTIGGRAPH_H 1

#include "DirectedGraph.h"
#include <istream>
#include <ostream>

class ContigGraph;

void readContigGraph(ContigGraph& graph, const std::string& path);

/** A contig graph is a directed graph with the property that
 * the edge (u,v) implies the existence of the edge (~v,~u).
 */
class ContigGraph : public DirectedGraph<> {
	typedef DirectedGraph<> DG;

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

	friend std::istream& operator>>(std::istream& in, ContigGraph& g);
	friend std::istream& operator<<(std::istream& in, ContigGraph& g);

  private:
	ContigGraph(const ContigGraph&);
};

/** Output a contig adjacency graph. */
static inline std::ostream& operator<<(std::ostream& out,
		const ContigGraph& g)
{
	for (ContigGraph::vertex_iterator v = g.begin();
			v != g.end(); ++v) {
		const ContigNode& id = g.vertex(*v);
		if (!id.sense())
			out << g.vertex(*v).id();
	   	out << "\t;";
		for (ContigGraph::out_edge_iterator e = v->begin();
				e != v->end(); ++e)
			out << ' ' << (g.target(*e) ^ id.sense());
		if (id.sense())
			out << '\n';
	}
	return out;
}

#endif
