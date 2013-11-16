#ifndef CONSTSTRING_H
#define CONSTSTRING_H 1

#include "Common/Hash.h"
#include "Common/HashFunction.h"
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

	bool operator==(const cstring& o) const
	{
		return m_p == o.m_p || strcmp(m_p, o.m_p) == 0;
	}

	bool operator<(const cstring& o) const
	{
		return m_p != o.m_p && strcmp(m_p, o.m_p) < 0;
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

#if __GXX_EXPERIMENTAL_CXX0X__
	const_string(const_string&& s) : cstring(s.m_p) { s.m_p = NULL; }
#endif

#if 0
	/* Should be like this, but... */
	const_string(const const_string& s)
		: cstring(strcpy(new char[s.size() + 1], s.c_str())) { }
#else
	/** Copy constructor.
	 * When a vector grows, libstdc++ calls the copy constructor for
	 * each element of the vector, which would invalidate any cstring
	 * that point to this const_string. To work around this issue, the
	 * new const_string gets the original data, and the old
	 * const_string gets the copy, which will probably be destructed
	 * soon. Making the copy is wasteful, but the C++ standard does
	 * not help us out here.
	 */
	const_string(const const_string& s) : cstring(s.c_str())
	{
		const_cast<const_string&>(s).m_p
			= strcpy(new char[s.size() + 1], s.c_str());
	}
#endif

	~const_string() { delete[] m_p; }

	const_string& operator=(const const_string& s)
	{
		assert(false);
		if (this == &s)
			return *this;
		assert(m_p != s.m_p);
		delete[] m_p;
		m_p = strcpy(new char[s.size() + 1], s.c_str());
		return *this;
	}

	void swap(const_string& s) { std::swap(m_p, s.m_p); }

  private:
	const_string();
	const_string(const char* s);
	const_string(const cstring&);
	bool operator==(const const_string& s);
};

namespace std {
	template <>
	inline void swap(const_string& a, const_string& b)
	{
		a.swap(b);
	}
}

NAMESPACE_STD_HASH_BEGIN
	template <> struct hash<cstring> {
		size_t operator()(const cstring& s) const
		{
			return hashmem(s.c_str(), strlen(s.c_str()));
		}
	};
NAMESPACE_STD_HASH_END

#endif
