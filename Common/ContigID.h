#ifndef CONTIGID_H
#define CONTIGID_H 1

#include "ConstString.h"
#include "Dictionary.h"
#include <cassert>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>

/** A contig ID is stored in memory as an integer, but is formatted as
 * a string using a static dictionary.
 */
class ContigID {
  public:
	ContigID() { }
	explicit ContigID(unsigned id) : m_id(id) { };
	explicit ContigID(const std::string& id)
		: m_id(s_dict.index(id)) { };

	/** Return the index of this ID. */
	operator unsigned() const { return m_id; }

	/** Return the string representation. */
	Dictionary::key_reference str() const
	{
		return s_dict.name(m_id);
	}

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

	friend std::ostream& operator <<(std::ostream& out,
			const ContigID& o)
	{
		return out << o.str();
	}

	static bool empty() { return s_dict.empty(); }
	static void lock() { s_dict.lock(); }
	static void unlock() { s_dict.unlock(); }
	static size_t count(const std::string& id)
	{
		return s_dict.count(id);
	}

	/** Create a new contig ID from s, which must be unique. */
	static ContigID insert(const std::string& s)
	{
		return ContigID(s_dict.insert(s));
	}

	/** Set the next contig ID returned by ContigID::create. */
	static void setNextContigID(cstring s)
	{
		std::istringstream iss((std::string)s);
		if (iss >> s_nextID) {
			assert(iss.eof());
			++s_nextID;
		} else
			s_nextID = 0;
	}

	/** Return a unique contig ID. */
	static ContigID create()
	{
		if (s_nextID == 0) {
			assert(!s_dict.empty());
			setNextContigID(s_dict.back());
		}
		std::ostringstream oss;
		oss << s_nextID++;
		return insert(oss.str());
	}

  private:
	/** The index. */
	unsigned m_id;

	/** The contig ID dictionary. */
	static Dictionary s_dict;

	/** The next unique contig ID. */
	static unsigned s_nextID;
};

#endif
