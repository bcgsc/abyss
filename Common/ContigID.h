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

/** Return the next unique contig name. */
static inline std::string createContigName()
{
	if (g_nextContigName == 0) {
		assert(!g_contigNames.empty());
		setNextContigName(g_contigNames.back());
	}
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
