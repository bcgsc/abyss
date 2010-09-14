#ifndef CONSTSTRING_H
#define CONSTSTRING_H 1

#include <cassert>
#include <cstring>
#include <ostream>
#include <string>

/** An immutable string that does not allocate resources. */
class cstring {
  public:
	cstring(const char* p) : m_p(p) { }
	cstring(const std::string& s) : m_p(s.c_str()) { }

	/** Return the size of this string. */
	size_t size() const { return strlen(m_p); }

	/** Return a null-terminated sequence of characters. */
	const char* c_str() const { return m_p; }

	operator const char*() const { return m_p; }

	bool operator<(const cstring& o) const
	{
		return strcmp(m_p, o.m_p) < 0;
	}

	friend std::ostream& operator<<(std::ostream& out,
			const cstring& o)
	{
		return out << o.m_p;
	}

  protected:
	const char* m_p;
};

/** An immutable string. */
class const_string : public cstring {
  public:
	const_string(const std::string& s)
		: cstring(strcpy(new char[s.size() + 1], s.c_str())) { }
	const_string(const const_string& s)
		: cstring(strcpy(new char[s.size() + 1], s.c_str())) { }

	~const_string() { delete[] m_p; }

	const_string& operator=(const const_string& s)
	{
		assert(false);
		delete[] m_p;
		m_p = strcpy(new char[s.size() + 1], s.c_str());
	}

	void swap(const_string& s) { std::swap(m_p, s.m_p); }

  private:
	const_string();
	const_string(const char* s);
	const_string(const cstring&);
};

namespace std {
	template <>
	inline void swap(const_string& a, const_string& b)
	{
		a.swap(b);
	}
}

#endif
