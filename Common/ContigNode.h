#ifndef CONTIGNODE_H
#define CONTIGNODE_H 1

#include "Dictionary.h"
#include "StringUtil.h"
#include <cassert>
#include <string>
#include <istream>
#include <ostream>

/** A dictionary of contig IDs. */
extern Dictionary g_contigIDs;

/** A tuple of a contig ID and an orientation. */
class ContigNode {
  public:
	ContigNode() { }
	ContigNode(unsigned id, bool sense) : m_node(id << 1 | sense) { };
	ContigNode(std::string id, bool sense)
		: m_node(g_contigIDs.serial(id) << 1 | sense) { };
	ContigNode(std::string id)
	{
		char c = chop(id);
		assert(c == '+' || c == '-');
		m_node = ContigNode(id, c == '-').m_node;
	}

	unsigned id() const { return m_node >> 1; }
	bool sense() const { return m_node & 1; }

	void flip() { m_node ^= 1; }

	bool operator ==(const ContigNode& o) const
	{
		return m_node == o.m_node;
	}

	bool operator <(const ContigNode& o) const
	{
		return m_node < o.m_node;
	}

	const ContigNode operator~() const
	{
		ContigNode o;
		o.m_node = this->m_node ^ 1;
		return o;
	}

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
		return out << g_contigIDs.key(o.id())
			<< (o.sense() ? '-' : '+');
	}

	// These functions are implemented in Overlap.
	unsigned outDegree() const;
	unsigned inDegree() const;
	const std::string sequence() const;

  private:
	unsigned m_node;
};

#endif
