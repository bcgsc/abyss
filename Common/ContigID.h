#ifndef CONTIGID_H
#define CONTIGID_H 1

#include "ConstString.h"
#include "Dictionary.h"
#include <boost/property_map/property_map.hpp>
#include <cassert>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>

/** The dictionary of contig names. */
extern Dictionary g_contigNames;

/** A contig ID is stored in memory as an integer, but is formatted as
 * a string using a static dictionary.
 */
class ContigID {
  public:
	ContigID() { }
	explicit ContigID(unsigned id) : m_id(id) { };
	explicit ContigID(const std::string& id)
		: m_id(g_contigNames.index(id)) { };

	/** Return the index of this ID. */
	operator unsigned() const { return m_id; }

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
		else
			o.m_id = 1<<29; // invalid value
		return in;
	}

	/** The insertion operator is not implemented to prevent
	 * unintentional use of the cast operator (unsigned).
	 */
	friend std::ostream& operator<<(std::ostream&, const ContigID&);

	/** Set the next contig ID returned by ContigID::create. */
	static void setNextContigID(cstring s)
	{
		std::istringstream iss((std::string)s);
		if (iss >> s_nextID && iss.eof())
			++s_nextID;
		else
			s_nextID = 0;
	}

	/** Return a unique contig ID. */
	static ContigID create()
	{
		if (s_nextID == 0) {
			assert(!g_contigNames.empty());
			setNextContigID(g_contigNames.back());
		}
		std::ostringstream oss;
		oss << s_nextID++;
		return ContigID(g_contigNames.insert(oss.str()));
	}

  private:
	/** The index. */
	unsigned m_id;

	/** The next unique contig ID. */
	static unsigned s_nextID;
};

/** A property map of a ContigID to an index. */
struct ContigIDIndexMap {
	typedef ContigID key_type;
	typedef unsigned value_type;
	typedef value_type reference;
	typedef boost::readable_property_map_tag category;
};

/** Return a numeric index of the specified contig. */
static inline
unsigned get(const ContigIDIndexMap&, ContigID u)
{
	return u;
}

#endif
