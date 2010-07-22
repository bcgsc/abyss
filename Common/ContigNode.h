#ifndef CONTIGNODE_H
#define CONTIGNODE_H 1

#include "config.h" // for WORDS_BIGENDIAN
#include "ContigID.h"
#include "StringUtil.h"
#include <cassert>
#include <cstdlib> // for strtoul
#include <string>
#include <istream>
#include <ostream>

/** A tuple of a contig ID and an orientation. */
class ContigNode {
  public:
	ContigNode() { }

#if WORDS_BIGENDIAN
	ContigNode(unsigned id, bool sense)
		: m_ambig(false), m_id(id), m_sense(sense) { }
	ContigNode(unsigned id, int sense)
		: m_ambig(false), m_id(id), m_sense(sense) { }
	ContigNode(std::string id, bool sense)
		: m_ambig(false),
		m_id(g_contigIDs.serial(id)), m_sense(sense) { }
#else
	ContigNode(unsigned id, bool sense)
		: m_sense(sense), m_id(id), m_ambig(false) { }
	ContigNode(unsigned id, int sense)
		: m_sense(sense), m_id(id), m_ambig(false) { }
	ContigNode(std::string id, bool sense)
		: m_sense(sense), m_id(g_contigIDs.serial(id)),
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
		assert(m_id < 100000);
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
	unsigned length() const;
	unsigned coverage() const;

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

#endif
