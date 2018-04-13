#ifndef CONTIGID_H
#define CONTIGID_H 1

#include "ConstString.h"
#include "Dictionary.h"
#include <cassert>
#include <ostream>
#include <sstream>
#include <string>

/** The dictionary of contig names. */
extern Dictionary g_contigNames;

/** The next unique contig name. */
extern unsigned g_nextContigName;

/** Set the next contig name returned by createContigName. */
static inline void setNextContigName(cstring s)
{
	std::istringstream iss((std::string)s);
	if (iss >> g_nextContigName && iss.eof())
		++g_nextContigName;
	else
		g_nextContigName = 0;
}

/**
 * Set the next contig name returned by createContigName
 * to one plus the current largest numeric contig name.
 */
static inline void setNextContigName()
{
	if (g_contigNames.empty()) {
		g_nextContigName = 0;
		return;
	}
	unsigned maxContigName = 0;
	for (unsigned i = 0; i < g_contigNames.size(); ++i) {
		cstring s = g_contigNames.getName(i);
		std::istringstream iss((std::string)s);
		unsigned contigName;
		if (iss >> contigName && iss.eof() && contigName > maxContigName)
			maxContigName = contigName;
	}
	g_nextContigName = 1 + maxContigName;
}

/** Return the next unique contig name. */
static inline std::string createContigName()
{
	if (g_nextContigName == 0)
		setNextContigName();
	std::ostringstream ss;
	ss << g_nextContigName++;
	return ss.str();
}

/** A contig index. */
class ContigID
{
  public:
	ContigID() { }
	explicit ContigID(unsigned index) : m_index(index) { };

	/** Return the index. */
	operator unsigned() const { return m_index; }

	/** The insertion operator is not implemented to prevent
	 * unintentional use of the cast operator (unsigned).
	 */
	friend std::ostream& operator<<(std::ostream&, const ContigID&);

  private:
	/** The index. */
	unsigned m_index;
};

#endif
