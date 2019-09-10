#ifndef CONTIGNODE_H
#define CONTIGNODE_H 1

#include "config.h" // for WORDS_BIGENDIAN
#include "ContigID.h"
#include "Graph/Properties.h"
#include "StringUtil.h"
#include <boost/property_map/property_map.hpp>
#include <cassert>
#include <cstdlib> // for strtoul
#include <iostream>
#include <stdint.h> // for intptr_t
#include <string>
#include <utility>

/** A vertex of a contig graph, which is a pair of a contig index and
 * an orientation.
 */
class ContigNode
{
public:

ContigNode() : m_index(0) { }

ContigNode(const ContigNode& o) : m_index(o.m_index) { }

/** Construct from a vertex index. */
explicit ContigNode(unsigned index) : m_index(index) { }

/** Construct from a contig index and an orientation. */
ContigNode(unsigned index, bool sense)
	: m_index(2 * index + sense) { }

/** Construct from a contig index and an orientation. */
ContigNode(unsigned index, int sense)
	: m_index(2 * index + sense)
{
	assert(sense == 0 || sense == 1);
}

/** Construct an ambiguous ContigNode. */
ContigNode(unsigned n, char c)
	: m_index(-(int)n)
{
	assert(n > 0);
	assert(c == 'N');
	(void)c;
}

bool operator==(const ContigNode& o) const
{
	return m_index == o.m_index;
}

bool operator!=(const ContigNode& o) const
{
	return m_index != o.m_index;
}

bool operator<(const ContigNode& o) const
{
	return m_index < o.m_index;
}

/** Return the complement of this vertex if sense is true. */
ContigNode operator^(bool sense) const
{
	assert(!ambiguous());
	return ContigNode(m_index ^ sense);
}

/** Copy constructors */
ContigNode(ContigNode&&) = default;
ContigNode& operator=(const ContigNode&) = default;
ContigNode& operator=(ContigNode&&) = default;

/** Return whether this ContigNode is ambiguous. */
bool ambiguous() const
{
	return m_index < 0;
}

/** Return the vertex index. */
unsigned index() const
{
	assert(!ambiguous());
	return m_index;
}

/** Return the contig index. */
unsigned id() const
{
	return ambiguous() ? m_index : m_index / 2;
}

/** Return the contig index as a ContigID. */
ContigID contigIndex() const
{
	assert(!ambiguous());
	return ContigID(id());
}

/** Return the orientation of this vertex. */
bool sense() const
{
	assert(!ambiguous());
	return m_index & 1;
}

/** Return the length in k-mer of this ambiguous contig. */
unsigned length() const
{
	assert(ambiguous());
	return -m_index;
}

/** Return the string of Ns. */
std::string ambiguousSequence() const
{
	assert(ambiguous());
	unsigned n = length();
	if (n > 100000) {
		std::cerr
			<< "warning: scaffold gap is longer than 100 kbp: "
			<< n << '\n';
	} else if (n > 1000000) {
		std::cerr << "error: scaffold gap is longer than 1 Mbp: "
			<< n << '\n';
		exit(EXIT_FAILURE);
	}
	return std::string(n, 'N');
}

/** Toggle the orientation of this vertex if sense is true. */
ContigNode& operator^=(bool sense)
{
	assert(!ambiguous());
	m_index ^= sense;
	return *this;
}

/** Increment this vertex index. */
ContigNode& operator++()
{
	assert(!ambiguous());
	++m_index;
	return *this;
}

private:
	int m_index;
};

/** Return the hash value of this ContigNode. */
static inline unsigned hash_value(const ContigNode& o)
{
	return o.index();
}

/** Vertex index property map of a ContigNode. */
struct ContigNodeIndexMap
	: boost::put_get_helper<unsigned, ContigNodeIndexMap>
{
	typedef ContigNode key_type;
	typedef unsigned value_type;
	typedef value_type reference;
	typedef boost::readable_property_map_tag category;

	reference operator[](const key_type& u) const
	{
		return u.index();
	}
};

/** Contig index property map of a ContigNode. */
struct ContigIndexMap
	: boost::put_get_helper<ContigID, ContigIndexMap>
{
	typedef ContigNode key_type;
	typedef ContigID value_type;
	typedef value_type reference;
	typedef boost::readable_property_map_tag category;

	reference operator[](const key_type& u) const
	{
		return u.contigIndex();
	}
};

/** Return the contig name of the specified vertex. */
template <typename Graph>
Dictionary::name_reference
get(vertex_contig_name_t, const Graph&, ContigNode u)
{
	assert(!u.ambiguous());
	return get(g_contigNames, u.id());
}

/** The string representation of a vertex name. */
struct VertexName : std::pair<intptr_t, char>
{
	typedef std::pair<intptr_t, char> Base;

	VertexName(const char* s, char c)
		: Base(reinterpret_cast<intptr_t>(s), c)
	{
		assert(c == '+' || c == '-');
	}

	VertexName(unsigned n, char c) : Base(n, c)
	{
		assert(c == 'N');
	}

	friend std::ostream& operator<<(
			std::ostream& out, const VertexName& o)
	{
		if (o.second == 'N')
			out << o.first;
		else
			out << reinterpret_cast<const char*>(o.first);
		return out << o.second;
	}
};

/** Return the name of the specified vertex. */
static inline VertexName get(const Dictionary& pmap, ContigNode u)
{
	return u.ambiguous()
		? VertexName(u.length(), 'N')
		: VertexName(get(pmap, u.id()), u.sense() ? '-' : '+');
}

/** Return the name of the specified vertex. */
template <typename Graph>
VertexName get(vertex_name_t, const Graph&, ContigNode u)
{
	return get(g_contigNames, u);
}

/** Set the name of the specified vertex. */
template <typename Graph>
void put(vertex_name_t, const Graph&, ContigNode u,
		const std::string& name)
{
	assert(!name.empty());
	char c = name[name.size() - 1];
	if (c == '+' || c == '-')
		put(g_contigNames, u.id(), name.substr(0, name.size() - 1));
	else
		put(g_contigNames, u.id(), name);
}

/** The string representation of an edge name. */
struct EdgeName : std::pair<VertexName, VertexName>
{
	typedef std::pair<VertexName, VertexName> Base;

	EdgeName(Base::first_type a, Base::second_type b)
		: Base(a, b) { }

	friend std::ostream& operator<<(
			std::ostream& out, const EdgeName& o)
	{
		return out << '"' << o.first << "\" -> \"" << o.second << '"';
	}
};

/** Return the name of the specified edge. */
static inline EdgeName get(const Dictionary& pmap,
		std::pair<ContigNode, ContigNode> e)
{
	return EdgeName(get(pmap, e.first), get(pmap, e.second));
}

/** Return the name of the specified edge. */
template <typename Graph>
EdgeName get(edge_name_t, const Graph& g,
		std::pair<ContigNode, ContigNode> e)
{
	return EdgeName(
			get(vertex_name, g, source(e, g)),
			get(vertex_name, g, target(e, g)));
}

/** Find the vertex with the specified name. */
static inline ContigNode find_vertex(
		std::string name, bool sense, const Dictionary& pmap)
{
	assert(!name.empty());
	return ContigNode(get(pmap, name), sense);
}

/** Find the vertex with the specified name. */
static inline ContigNode find_vertex(
		std::string name, const Dictionary& pmap)
{
	assert(!name.empty());
	char c = chop(name);
	assert(c == '+' || c == '-' || c == 'N');
	return c == 'N'
		? ContigNode(strtoul(name.c_str(), NULL, 0), 'N')
		: find_vertex(name, c == '-', pmap);
}

#endif
