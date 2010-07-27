#ifndef CONTIGID_H
#define CONTIGID_H 1

#include "Dictionary.h"
#include <istream>
#include <ostream>
#include <string>

extern Dictionary g_contigIDs;

/** A contig ID is represented by a numeric serial number, but is
 * formatted as a string.
 */
class ContigID {
  public:
	ContigID() : m_id(0) { }
	explicit ContigID(unsigned id) : m_id(id) { };
	explicit ContigID(std::string id)
		: m_id(g_contigIDs.serial(id)) { };

	/** Return the numeric serial number. */
	operator unsigned() const { return m_id; }

	/** Return the string representation. */
	const std::string& str() const { return g_contigIDs.key(m_id); }

	bool operator ==(const ContigID& o) const
	{
		return m_id == o.m_id;
	}

	bool operator <(const ContigID& o) const
	{
		return m_id < o.m_id;
	}

	friend std::istream& operator >>(std::istream& in,
			ContigID& o)
	{
		std::string s;
		if (in >> s)
			o = ContigID(s);
		return in;
	}

	friend std::ostream& operator <<(std::ostream& out,
			const ContigID& o)
	{
		return out << o.str();
	}

  private:
	/** The numeric serial number. */
	unsigned m_id;
};

#endif
