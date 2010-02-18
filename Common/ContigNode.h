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
	ContigNode(unsigned id, bool sense) : m_id(id), m_sense(sense) { };
	ContigNode(std::string id, bool sense)
		: m_id(g_contigIDs.serial(id)), m_sense(sense) { };
	ContigNode(std::string id)
	{
		char c = chop(id);
		assert(c == '+' || c == '-');
		*this = ContigNode(id, c == '-');
	}

	unsigned id() const { return m_id; }
	bool sense() const { return m_sense; }

	void flip() { m_sense = !m_sense; }

	bool operator ==(const ContigNode& o) const
	{
		return m_id == o.m_id && m_sense == o.m_sense;
	}

	bool operator <(const ContigNode& o) const
	{
		return hash() < o.hash();
	}

	const ContigNode operator~() const
	{
		return ContigNode(m_id, !m_sense);
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
	unsigned length() const;
	const std::string sequence() const;

  private:
	unsigned hash() const { return m_id << 1 | m_sense; }
	unsigned m_id:31;
	unsigned m_sense:1;
};

#endif
