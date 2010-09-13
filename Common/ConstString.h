#ifndef CONSTSTRING_H
#define CONSTSTRING_H 1

#include <cassert>
#include <cstring>
#include <ostream>
#include <string>

/** An immutable string. */
class const_string
{
  public:
	const_string(const char* p) : m_p(p) { }
	const_string(const std::string& s) : m_p(s.c_str()) { }

	/** Make a copy of this string. Use free to release it. */
	const_string clone() const
	{
		return strcpy(new char[strlen(m_p) + 1], m_p);
	}

	/** Release the resources allocated using clone. */
	void free()
	{
		assert(m_p != NULL);
		delete[] m_p;
		m_p = NULL;
	}

	/** Return a null-terminated sequence of characters. */
	const char* c_str() const { return m_p; }

	operator const char*() const { return m_p; }

	bool operator<(const const_string& o) const
	{
		return strcmp(m_p, o.m_p) < 0;
	}

	friend std::ostream& operator<<(std::ostream& out,
			const const_string& o)
	{
		return out << o.m_p;
	}

  private:
	const char* m_p;
};

#endif
