#ifndef DOTIO_H
#define DOTIO_H 1

#include "ContigID.h" // for g_contigNames.lock
#include "Graph/Options.h"
#include "IOUtil.h"
#include <boost/graph/graph_traits.hpp>
#include <cassert>
#include <cstdlib> // for exit
#include <istream>
#include <ostream>

using boost::graph_traits;

template <typename Graph, typename VertexProp>
void write_vertex(std::ostream& out, const Graph& g,
		typename graph_traits<Graph>::vertex_descriptor u,
		const VertexProp*)
{
	out << '"' << get(vertex_name, g, u) << "\""
		" [" << get(vertex_bundle, g, u) << "]\n";
}

template <typename Graph>
void write_vertex(std::ostream&, const Graph&,
		typename graph_traits<Graph>::vertex_descriptor,
		const no_property*)
{
}

template <typename Graph, typename EdgeProp>
void write_edges(std::ostream& out, const Graph& g,
		typename graph_traits<Graph>::vertex_descriptor u,
		const EdgeProp*)
{
	typedef typename graph_traits<Graph>::out_edge_iterator
		out_edge_iterator;
	typedef typename edge_property<Graph>::type edge_property_type;
	std::pair<out_edge_iterator, out_edge_iterator>
		adj = out_edges(u, g);
	for (out_edge_iterator e = adj.first; e != adj.second; ++e) {
		assert(!get(vertex_removed, g, target(*e, g)));
		out << get(edge_name, g, *e);
		const edge_property_type& ep = get(edge_bundle, g, e);
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

/** Output a GraphViz dot graph. */
template <typename Graph>
std::ostream& write_dot(std::ostream& out, const Graph& g)
{
	return write_dot(out, g, "adj");
}

/** Output a GraphViz dot graph. */
template <typename Graph>
std::ostream& write_dot(std::ostream& out, const Graph& g,
	const std::string& graphName)
{
	typedef typename graph_traits<Graph>::vertex_iterator
		vertex_iterator;
	typedef typename vertex_property<Graph>::type
		vertex_property_type;
	typedef typename edge_property<Graph>::type edge_property_type;

	out << "digraph " << graphName << " {\n";

	// Graph properties
	if (opt::k > 0)
		out << "graph [k=" << opt::k << "]\n"
			"edge [d=" << -int(opt::k - 1) << "]\n";

	// Vertices
	std::pair<vertex_iterator, vertex_iterator> vit = vertices(g);
	for (vertex_iterator u = vit.first; u != vit.second; ++u) {
		if (get(vertex_removed, g, *u))
			continue;
		write_vertex(out, g, *u, (vertex_property_type*)NULL);
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

/** legal characters for an unquoted node ID */
static const uint8_t GRAPHVIZ_ID_CHARS[256] = {
	0, 0, 0, 0, 0, 0, 0, 0, // 000-007:
	0, 0, 0, 0, 0, 0, 0, 0, // 008-015:
	0, 0, 0, 0, 0, 0, 0, 0, // 016-023:
	0, 0, 0, 0, 0, 0, 0, 0, // 024-031:
	0, 0, 0, 0, 0, 0, 0, 0, // 032-039:    !, ", #, $, %, &, ',
	0, 0, 0, 1, 0, 1, 1, 0, // 040-047: (, ), *, +, ,, -, ., /,
	1, 1, 1, 1, 1, 1, 1, 1, // 048-055: 0, 1, 2, 3, 4, 5, 6, 7,
	1, 1, 0, 0, 0, 0, 0, 0, // 056-063: 8, 9, :, ;, <, =, >, ?,
	0, 1, 1, 1, 1, 1, 1, 1, // 064-071: @, A, B, C, D, E, F, G,
	1, 1, 1, 1, 1, 1, 1, 1, // 072-079: H, I, J, K, L, M, N, O,
	1, 1, 1, 1, 1, 1, 1, 1, // 080-087: P, Q, R, S, T, U, V, W,
	1, 1, 1, 0, 0, 0, 0, 1, // 088-095: X, Y, Z, [, \, ], ^, _,
	0, 1, 1, 1, 1, 1, 1, 1, // 096-103: `, a, b, c, d, e, f, g,
	1, 1, 1, 1, 1, 1, 1, 1, // 104-111: h, i, j, k, l, m, n, o,
	1, 1, 1, 1, 1, 1, 1, 1, // 112-119: p, q, r, s, t, u, v, w,
	1, 1, 1, 0, 0, 0, 0, 0, // 120-127: x, y, z, {, |, }, ~,
	0, 0, 0, 0, 0, 0, 0, 0, // 128-135:
	0, 0, 0, 0, 0, 0, 0, 0, // 136-143:
	0, 0, 0, 0, 0, 0, 0, 0, // 144-151:
	0, 0, 0, 0, 0, 0, 0, 0, // 152-159:
	0, 0, 0, 0, 0, 0, 0, 0, // 160-167:
	0, 0, 0, 0, 0, 0, 0, 0, // 168-175:
	0, 0, 0, 0, 0, 0, 0, 0, // 176-183:
	0, 0, 0, 0, 0, 0, 0, 0, // 184-191:
	0, 0, 0, 0, 0, 0, 0, 0, // 192-199:
	0, 0, 0, 0, 0, 0, 0, 0, // 200-207:
	0, 0, 0, 0, 0, 0, 0, 0, // 208-215:
	0, 0, 0, 0, 0, 0, 0, 0, // 216-223:
	0, 0, 0, 0, 0, 0, 0, 0, // 224-231:
	0, 0, 0, 0, 0, 0, 0, 0, // 232-239:
	0, 0, 0, 0, 0, 0, 0, 0, // 240-247:
	0, 0, 0, 0, 0, 0, 0, 0  // 248-255:
};

/** Read a string delimited by double quotes. */
static inline
std::istream& read_dot_name(std::istream& in, std::string& s)
{
	in >> std::ws;
	if (in.peek() == '"') {
		// quoted node ID
		in.get();
		getline(in, s, '"');
	} else if (GRAPHVIZ_ID_CHARS[in.peek()]) {
		// unquoted node ID
		std::ostringstream ss;
		while (GRAPHVIZ_ID_CHARS[in.peek()])
			ss << (char)in.get();
		assert(ss.good());
		s = ss.str();
	} else {
		in.clear(std::ios::failbit);
	}
	return in;
}

/** Read a vertex descriptor delimited by double quotes. */
template <typename Graph>
std::istream& read_dot_id(std::istream& in, const Graph& g,
		typename graph_traits<Graph>::vertex_descriptor& u)
{
	std::string s;
	if (read_dot_name(in, s))
		u = find_vertex(s, g);
	return in;
}

/** Read a GraphViz dot graph.
 * @param betterEP handle parallel edges
 */
template <typename Graph, typename BetterEP>
std::istream& read_dot(std::istream& in, Graph& g, BetterEP betterEP)
{
	assert(in);
	typedef typename graph_traits<Graph>::vertex_descriptor
		vertex_descriptor;
	typedef typename vertex_property<Graph>::type
		vertex_property_type;
	typedef typename graph_traits<Graph>::edge_descriptor
		edge_descriptor;
	typedef typename edge_property<Graph>::type edge_property_type;
	typedef typename graph_traits<Graph>::directed_category
		directed_category;

	bool isDirectedGraph = boost::detail::is_directed(directed_category());

	// Add vertices if this graph is empty.
	bool addVertices = num_vertices(g) == 0;

	// Graph properties
	in >> expect(isDirectedGraph ? "digraph" : "graph");
	in >> Ignore('{');
	assert(in);

	edge_property_type defaultEdgeProp;
	for (bool done = false; !done && in >> std::ws;) {
		switch (in.peek()) {
		  case 'g': {
			// Graph Properties
			in >> expect("graph [ ");
			if (in.peek() == 'k') {
				unsigned k;
				in >> expect("k =") >> k;
				assert(in);
				if (opt::k > 0)
					assert(k == opt::k);
				opt::k = k;
			}
			in >> Ignore(']');
			break;
		  }
		  case 'e': // edge
			// Default edge properties
			in >> expect("edge [") >> defaultEdgeProp >> Ignore(']');
			assert(in);
			break;
		  default:
			done = true;
			break;
		}
		if (in >> std::ws && in.peek() == ';')
			in.get();
	}

	for (std::string uname; read_dot_name(in, uname);) {
		char c;
		in >> c;
		assert(in);
		if (c == ';') {
			// Vertex
			if (addVertices) {
				vertex_descriptor u = add_vertex(g);
				put(vertex_name, g, u, uname);
			} else {
				vertex_descriptor u = find_vertex(uname, g);
				assert(get(vertex_index, g, u) < num_vertices(g));
				(void)u;
			}
		} else if (c == '[') {
			// Vertex properties
			vertex_property_type vp;
			in >> vp >> Ignore(']');
			assert(in);
			if (addVertices) {
				vertex_descriptor u = add_vertex(vp, g);
				put(vertex_name, g, u, uname);
			} else {
				vertex_descriptor u = find_vertex(uname, g);
				assert(get(vertex_index, g, u) < num_vertices(g));
				if (g[u] != vp) {
					std::cerr << "error: "
						"vertex properties do not agree: "
						"\"" << uname << "\" "
						"[" << g[u] << "] [" << vp << "]\n";
					exit(EXIT_FAILURE);
				}
			}
		} else if (c == '-') {
			// Edge
			in >> expect(isDirectedGraph ? ">" : "-");
			assert(in);
			g_contigNames.lock();

			vertex_descriptor u = find_vertex(uname, g);
			if (in >> std::ws && in.peek() == '{') {
				// Subgraph
				in >> expect("{");
				for (vertex_descriptor v; read_dot_id(in, g, v);)
					add_edge(u, v, defaultEdgeProp, g);
				if (in.fail())
					in.clear();
				in >> expect(" }");
				assert(in);
			} else {
				vertex_descriptor v;
				if (!read_dot_id(in, g, v)) {
					in.clear();
					std::cerr << "error: Expected `\"' and saw `"
						<< (char)in.peek() << "'.\n";
					exit(EXIT_FAILURE);
				}
				assert(in);

				edge_property_type ep = defaultEdgeProp;
				if (in >> std::ws && in.peek() == '[') {
					// Edge properties
					in >> expect("[") >> ep >> Ignore(']');
					assert(in);
				}

				edge_descriptor e;
				bool found;
				boost::tie(e, found) = edge(u, v, g);
				if (found) {
					// Parallel edge
					edge_property_type& ref = g[e];
					ref = betterEP(ref, ep);
				} else
					add_edge(u, v, ep, g);
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
