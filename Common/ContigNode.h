#ifndef CONTIGNODE_H
#define CONTIGNODE_H 1

#include "config.h" // for WORDS_BIGENDIAN
#include "ContigID.h"
#include "Graph/Properties.h" // for vertex_index_t
#include "StringUtil.h"
#include <boost/property_map/property_map.hpp>
#include <cassert>
#include <cstdlib> // for strtoul
#include <iostream>
#include <string>

/** A tuple of a contig ID and an orientation. */
class ContigNode {
  public:
	ContigNode() { }

#if WORDS_BIGENDIAN
	ContigNode(unsigned id, bool sense)
		: m_ambig(false), m_id(id), m_sense(sense) { }
	ContigNode(unsigned id, int sense)
		: m_ambig(false), m_id(id), m_sense(sense) { }
	ContigNode(const std::string& id, bool sense)
		: m_ambig(false),
		m_id(ContigID(id)), m_sense(sense) { }
#else
	ContigNode(unsigned id, bool sense)
		: m_sense(sense), m_id(id), m_ambig(false) { }
	ContigNode(unsigned id, int sense)
		: m_sense(sense), m_id(id), m_ambig(false) { }
	ContigNode(const std::string& id, bool sense)
		: m_sense(sense), m_id(ContigID(id)),
		m_ambig(false) { }
#endif

	explicit ContigNode(unsigned i) : m_int(i) { }

	/** Create an ambiguous contig. */
	ContigNode(unsigned n, char c)
#if WORDS_BIGENDIAN
		: m_ambig(true), m_id(n), m_sense(false)
#else
		: m_sense(false), m_id(n), m_ambig(true)
#endif
	{
		assert(c == 'N');
		(void)c;
		assert(n > 0);
	}

	ContigNode(std::string id)
	{
		char c = chop(id);
		assert(c == '+' || c == '-' || c == 'N');
		*this = c == 'N'
			? ContigNode(strtoul(id.c_str(), NULL, 0), 'N')
			: ContigNode(id, c == '-');
	}

	bool ambiguous() const { return m_ambig; }
	unsigned id() const { return m_ambig ? -m_id : m_id; }
	bool sense() const { assert(!m_ambig); return m_sense; }

	std::string ambiguousSequence() const
	{
		assert(m_ambig);
		if (m_id > 100000) {
			std::cerr
				<< "warning: scaffold gap is longer than 100 kbp: "
				<< *this << '\n';
		} else if (m_id > 1000000) {
			std::cerr << "error: scaffold gap is longer than 1 Mbp: "
				<< *this << '\n';
			exit(EXIT_FAILURE);
		}
		return std::string(m_id, 'N');
	}

	void flip() { if (!m_ambig) m_sense = !m_sense; }

	/** Return the contig ID. */
	operator ContigID() const
	{
		assert(!m_ambig);
		return ContigID(m_id);
	}

	bool operator ==(const ContigNode& o) const
	{
		return hash() == o.hash();
	}

	bool operator !=(const ContigNode& o) const
	{
		return !(*this == o);
	}

	bool operator <(const ContigNode& o) const
	{
		return hash() < o.hash();
	}

	const ContigNode operator~() const
	{
		assert(!m_ambig);
		return ContigNode(m_id, !m_sense);
	}

	const ContigNode operator^(bool flip) const
	{
		return flip ? ~*this : *this;
	}

	ContigNode& operator++() { ++m_int; return *this; }

	friend std::istream& operator >>(std::istream& in,
			ContigNode& o)
	{
		std::string s;
		if (in >> s)
			o = ContigNode(s);
		return in;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const ContigNode& o)
	{
		if (o.ambiguous())
			return out << o.m_id << 'N';
		else
			return out << ContigID(o.id()) << (o.sense() ? '-' : '+');
	}

	/** Return the length of this ambiguous contig in k-mer. */
	unsigned length() const { assert(m_ambig); return m_id; }

	/** Return a value that can be used as an index of an array. */
	unsigned index() const
	{
		assert(!m_ambig);
		return hash();
	}

  private:
	/** Return the hash value of this contig. */
	unsigned hash() const { return m_int; }

	union {
		unsigned m_int;
		struct {
#if WORDS_BIGENDIAN
			unsigned m_ambig:1;
			unsigned m_id:30;
			unsigned m_sense:1;
#else
			unsigned m_sense:1;
			unsigned m_id:30;
			unsigned m_ambig:1;
#endif
		};
	};
};

/** Return the hash value of this ContigNode. */
static inline unsigned hash_value(const ContigNode& o)
{
	return o.index();
}

/** Return a numeric index of the specified vertex. */
static inline unsigned index(const ContigNode& o)
{
	return o.index();
}

/** Vertex index property map of a ContigNode. */
struct ContigNodeIndexMap {
	typedef ContigNode key_type;
	typedef unsigned value_type;
	typedef value_type reference;
	typedef boost::readable_property_map_tag category;

	reference operator[](const key_type& u) const
	{
		return u.index();
	}
};

/** Return a numeric index of the specified vertex. */
static inline
unsigned get(ContigNodeIndexMap, ContigNode u)
{
	return u.index();
}

/** Return a numeric index of the specified vertex. */
template <typename Graph>
unsigned get(vertex_index_t, const Graph&, ContigNode u)
{
	return u.index();
}

/** Return the sense of the specified vertex. */
template <typename Graph>
bool get(vertex_sense_t, const Graph&, ContigNode u)
{
	return u.sense();
}

/** Set the name of the specified vertex. */
template <typename Graph>
void put(vertex_name_t, const Graph&, ContigNode u,
		const std::string& name)
{
	assert(!name.empty());
	char c = name[name.size() - 1];
	if (c == '+' || c == '-')
		ContigID::put(u, std::string(name, 0, name.size() - 1));
	else
		ContigID::put(u, name);
}

#endif
