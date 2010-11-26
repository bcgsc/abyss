#ifndef DOTIO_H
#define DOTIO_H 1

#include "Graph/Options.h"
#include "Graph.h"
#include "IOUtil.h"
#include <cassert>
#include <cstdlib> // for exit
#include <istream>
#include <ostream>

template <typename V, typename VertexProp>
void write_vertex(std::ostream& out,
		V u, const VertexProp& vp)
{
	out << '"' << u << "\" [" << vp << "]\n";
}

template <typename V>
void write_vertex(std::ostream&, V, no_property)
{
}

template <typename Graph, typename EdgeProp>
void write_edges(std::ostream& out, const Graph& g,
		typename graph_traits<Graph>::vertex_descriptor u,
		const EdgeProp*)
{
	typedef typename graph_traits<Graph>::adjacency_iterator
		adjacency_iterator;
	typedef typename graph_traits<Graph>::edge_descriptor
		edge_descriptor;
	typedef typename edge_property<Graph>::type edge_property_type;
	std::pair<adjacency_iterator, adjacency_iterator>
		adj = adjacent_vertices(u, g);
	for (adjacency_iterator v = adj.first; v != adj.second; ++v) {
		assert(!get(vertex_removed, g, *v));
		out << '"' << u << "\" -> \"" << *v << '"';
		const edge_property_type& ep = get(edge_bundle, g,
				edge_descriptor(u, *v));
		if (!(ep == edge_property_type()))
			out << " [" << ep << ']';
		out << '\n';
	}
}

template <typename Graph>
void write_edges(std::ostream& out, const Graph& g,
		typename graph_traits<Graph>::vertex_descriptor u,
		const no_property*)
{
	typedef typename graph_traits<Graph>::adjacency_iterator
		adjacency_iterator;
	unsigned outdeg = out_degree(u, g);
	if (outdeg == 0)
		return;
	out << '"' << u << "\" ->";
	if (outdeg > 1)
		out << " {";
	std::pair<adjacency_iterator, adjacency_iterator>
		adj = adjacent_vertices(u, g);
	for (adjacency_iterator v = adj.first; v != adj.second; ++v)
		out << " \"" << *v << '"';
	if (outdeg > 1)
		out << " }";
	out << '\n';
}

/** Output a GraphViz dot graph. */
template <typename Graph>
std::ostream& write_dot(std::ostream& out, const Graph& g)
{
	typedef typename graph_traits<Graph>::vertex_iterator
		vertex_iterator;
	typedef typename edge_property<Graph>::type edge_property_type;

	out << "digraph adj {\n";

	// Graph properties
	if (opt::k > 0)
		out << "graph [k=" << opt::k << "]\n"
			"edge [d=" << -int(opt::k - 1) << "]\n";

	// Vertices
	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u) {
		if (get(vertex_removed, g, *u))
			continue;
		write_vertex(out, *u, get(vertex_bundle, g, *u));
	}

	// Edges
	for (vertex_iterator u = vit.first; u != vit.second; ++u) {
		if (get(vertex_removed, g, *u))
			continue;
		write_edges(out, g, *u, (edge_property_type*)NULL);
	}

	return out << "}\n";
}

/** Output a GraphViz dot graph. */
template <typename Graph>
struct DotWriter
{
	const Graph& g;
	DotWriter(const Graph& g) : g(g) { }
	friend std::ostream& operator<<(std::ostream& out,
			const DotWriter& o)
	{
		return write_dot<Graph>(out, o.g);
	}
};

/** Output a GraphViz dot graph. */
template <typename Graph>
DotWriter<Graph> dot_writer(const Graph& g)
{
	return DotWriter<Graph>(g);
}

/** Read an id delimited by double quotes. */
template <typename VertexDescriptor>
std::istream& read_dot_id(std::istream& in, VertexDescriptor& u)
{
	if (in >> std::ws && in.peek() == '"') {
		in.get();
		std::string s;
		getline(in, s, '"');
		std::istringstream ss(s);
		ss >> u;
		assert(ss.eof());
	} else
		in.clear(std::ios::failbit);
	return in;
}

template <typename Graph>
std::istream& read_dot(std::istream& in, Graph& g)
{
	assert(in);
	assert(num_vertices(g) == 0);

	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	typedef typename vertex_property<Graph>::type
		vertex_property_type;
	typedef typename edge_property<Graph>::type edge_property_type;

	// Graph properties
	in >> expect("digraph") >> ignore('{');
	assert(in);

	edge_property_type defaultEdgeProp;
	for (bool done = false; !done && in >> std::ws;) {
		switch (in.peek()) {
		  case 'g': {
			// Graph Properties
			unsigned k;
			in >> expect("graph [ k =") >> k >> ignore(']');
			assert(in);
			if (opt::k > 0)
				assert(k == opt::k);
			opt::k = k;
			break;
		  }
		  case 'e': // edge
			// Default edge properties
			in >> expect("edge [") >> defaultEdgeProp >> ignore(']');
			assert(in);
			break;
		  default:
			done = true;
			break;
		}
		if (in >> std::ws && in.peek() == ';')
			in.get();
	}

	vertex_descriptor u;
	while (read_dot_id(in, u)) {
		char c;
		in >> c;
		assert(in);
		if (c == ';') {
			// Vertex
			vertex_descriptor x = add_vertex(g);
			assert(x == u);
		} else if (c == '[') {
			// Vertex properties
			vertex_property_type vp;
			in >> vp >> ignore(']');
			assert(in);
			vertex_descriptor x = add_vertex(vp, g);
			assert(x == u);
		} else if (c == '-') {
			// Edge
			in >> expect(">");
			assert(in);
			ContigID::lock();

			if (in >> std::ws && in.peek() == '{') {
				// Subgraph
				in >> expect("{");
				for (vertex_descriptor v; read_dot_id(in, v);)
					add_edge(u, v, defaultEdgeProp, g);
				if (in.fail())
					in.clear();
				in >> expect(" }");
				assert(in);
			} else {
				vertex_descriptor v;
				read_dot_id(in, v);
				if (in.fail()) {
					in.clear();
					std::cerr << "error: Expected `\"' and saw `"
						<< (char)in.peek() << "'.\n";
					exit(EXIT_FAILURE);
				}
				assert(in);
				if (in >> std::ws && in.peek() == '[') {
					// Edge properties
					edge_property_type ep;
					in >> expect("[") >> ep >> ignore(']');
					assert(in);
					add_edge(u, v, ep, g);
				} else
					add_edge(u, v, defaultEdgeProp, g);
			}
		} else {
			std::cerr << "error: Expected `[' or `->' and saw `"
				<< c << "'.\n";
			exit(EXIT_FAILURE);
		}
		if (in >> std::ws && in.peek() == ';')
			in.get();
	}
	if (in.fail())
		in.clear();

	// Check for the closing brace.
	in >> expect("}") >> std::ws;
	assert(in.eof());

	assert(num_vertices(g) > 0);
	return in;
}

#endif
