#ifndef CONTIGNODE_H
#define CONTIGNODE_H 1

#include "Dictionary.h"
#include "StringUtil.h"
#include <cassert>
#include <cstdlib> // for strtoul
#include <string>
#include <istream>
#include <ostream>

/** A dictionary of contig IDs. */
extern Dictionary g_contigIDs;

/** A tuple of a contig ID and an orientation. */
class ContigNode {
  public:
	ContigNode() { }
	ContigNode(unsigned id, bool sense)
		: m_ambig(false), m_id(id), m_sense(sense) { }
	ContigNode(std::string id, bool sense)
		: m_ambig(false),
		m_id(g_contigIDs.serial(id)), m_sense(sense) { }

	/** Create an ambiguous contig. */
	explicit ContigNode(unsigned n)
		: m_ambig(true), m_id(n), m_sense(false) { assert(n > 0); }

	ContigNode(std::string id)
	{
		char c = chop(id);
		assert(c == '+' || c == '-' || c == 'N');
		*this = c == 'N' ? ContigNode(strtoul(id.c_str(), NULL, 0))
			: ContigNode(id, c == '-');
	}

	bool ambiguous() const { return m_ambig; }
	unsigned id() const { return m_ambig ? -m_id : m_id; }
	bool sense() const { assert(!m_ambig); return m_sense; }

	std::string ambiguousSequence() const
	{
		assert(m_ambig);
		assert(m_id < 100000);
		return std::string(m_id, 'N');
	}

	void flip() { if (!m_ambig) m_sense = !m_sense; }

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
			return out << g_contigIDs.key(o.id())
				<< (o.sense() ? '-' : '+');
	}

	// These functions are implemented elsewhere.
	unsigned outDegree() const;
	unsigned inDegree() const;
	unsigned length() const;
	const std::string sequence() const;
	unsigned coverage() const;

  private:
	unsigned hash() const
	{
		return m_ambig << 31 | m_id << 1 | m_sense;
	}

	unsigned m_ambig:1;
	unsigned m_id:30;
	unsigned m_sense:1;
};

#endif
