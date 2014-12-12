#ifndef ASSEMBLY_DOTWRITER_H
#define ASSEMBLY_DOTWRITER_H 1

/** Written by Shaun Jackman <sjackman@bcgsc.ca>. */

#include "Common/UnorderedMap.h"
#include "Graph/ContigGraphAlgorithms.h"
#include <cassert>
#include <iostream>
#include <sstream>

class DotWriter
{
private:
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;
	typedef graph_traits<Graph>::vertex_iterator Vit;
	typedef graph_traits<Graph>::adjacency_iterator Ait;
	typedef std::string VertexName;

	DotWriter() : m_id(0) { }

/** Complement the specified name. */
static std::string
complementName(std::string s)
{
	char c = *s.rbegin();
	assert(c == '+' || c == '-');
	*s.rbegin() = c == '+' ? '-' : '+';
	return s;
}

/** Return the name of the specified vertex. */
const VertexName&
getName(const V& u) const
{
	Names::const_iterator it = m_names.find(u);
	if (it == m_names.end())
		std::cerr << "error: cannot find vertex " << u << '\n';
	assert(it != m_names.end());
	return it->second;
}

/** Set the name of the specified vertex. */
void setName(const V& u, const VertexName& uname)
{
	std::pair<Names::const_iterator, bool> inserted
		= m_names.insert(Names::value_type(u, uname));
	if (!inserted.second)
		std::cerr << "error: duplicate vertex " << u << '\n';
	assert(inserted.second);
	(void)inserted;
}

/** Return whether this vertex is contiguous and not palindromic. */
bool contiguousOut(const Graph& g, const V& u)
{
	return contiguous_out(g, u)
		&& !u.isPalindrome()
		&& !u.isPalindrome(SENSE)
		&& !(*adjacent_vertices(u, g).first).isPalindrome();
}

/** Return the in-adjacent vertex. */
V adjacentVertexIn(const V& u, const Graph& g)
{
	return source(*in_edges(u, g).first, g);
}

/** Return whether this vertex is contiguous and not palindromic. */
bool contiguousIn(const Graph& g, const V& u)
{
	return contiguous_in(g, u)
		&& !u.isPalindrome()
		&& !u.isPalindrome(ANTISENSE)
		&& !adjacentVertexIn(u, g).isPalindrome();
}

/** Write out the specified contig. */
void writeContig(std::ostream& out, const Graph& g, const V& u)
{
	unsigned n = 0;
	unsigned c = 0;
	V v;
	for (v = u; contiguousOut(g, v); v = *adjacent_vertices(v, g).first) {
		++n;
		c += get(vertex_coverage, g, v);
	}
	++n;
	c += get(vertex_coverage, g, v);

	// Output the canonical orientation of the contig.
	V vrc = get(vertex_complement, g, v);
	if (vrc < u)
		return;

	std::ostringstream ss;
	ss << m_id << '+';
	VertexName uname = ss.str();
	VertexName vname(uname);
	*vname.rbegin() = '-';
	++m_id;

	// Reorient the contig to agree with assembleContig.
	bool rc;
	Graph::const_iterator it = g.find(u, rc);
	assert(it != g.end());
	(void)it;
	if (rc)
		std::swap(uname, vname);

	setName(u, uname);
	if (u == vrc) {
		// Palindrome
		assert(n == 1);
	} else
		setName(vrc, vname);

	unsigned l = n + V::length() - 1;
	out << '"' << uname << "\" [l=" << l << " C=" << c << "]\n"
		"\"" << vname << "\" [l=" << l << " C=" << c << "]\n";
}

/** Write out the contigs that split at the specified sequence. */
void
writeEdges(std::ostream& out, const Graph& g,
		const V& u, const VertexName& uname) const
{
	if (out_degree(u, g) == 0)
		return;
	out << '"' << uname << "\" -> {";
	std::pair<Ait, Ait> adj = adjacent_vertices(u, g);
	for (Ait vit = adj.first; vit != adj.second; ++vit) {
		V v = *vit;
		const VertexName& vname = getName(v);
		out << " \"" << vname << '"';
		if (v.isPalindrome())
			out << " \"" << complementName(vname) << '"';
	}
	out << " }\n";
}

/** Output the edges of the specified vertex. */
void
writeEdges(std::ostream& out, const Graph& g, const V& u) const
{
	std::string uname = complementName(getName(get(vertex_complement, g, u)));
	writeEdges(out, g, u, uname);
	if (u.isPalindrome()) {
		uname = complementName(uname);
		writeEdges(out, g, u, uname);
	}
}

/** Write out a dot graph for the specified collection. */
void writeGraph(std::ostream& out, const Graph& g)
{
	out << "digraph g {\n"
		"graph [k=" << V::length() << "]\n"
		"edge [d=" << -int(V::length() - 1) << "]\n";
	std::pair<Vit, Vit> uits = vertices(g);

	// Output the vertices.
	for (Vit uit = uits.first; uit != uits.second; ++uit) {
		V u = *uit;
		if (get(vertex_removed, g, u))
			continue;
		if (!contiguousIn(g, u))
			writeContig(out, g, u);
		// Skip the second occurence of the palindrome.
		if (u.isPalindrome()) {
			assert(uit != uits.second);
			++uit;
		}
	}

	// Output the edges.
	for (Vit uit = uits.first; uit != uits.second; ++uit) {
		V u = *uit;
		if (get(vertex_removed, g, u))
			continue;
		if (!contiguousOut(g, u))
			writeEdges(out, g, u);
		// Skip the second occurence of the palindrome.
		if (u.isPalindrome()) {
			assert(uit != uits.second);
			++uit;
		}
	}

	out << "}" << std::endl;
}

public:

/** Write out a dot graph for the specified collection. */
static
void write(std::ostream& out, const Graph& g)
{
	DotWriter dotWriter;
	dotWriter.writeGraph(out, g);
}

private:
	typedef unordered_map<V, VertexName> Names;

	/** A map of terminal k-mers to contig names. */
	Names m_names;

	/** The current contig name. */
	unsigned m_id;
};

#endif
