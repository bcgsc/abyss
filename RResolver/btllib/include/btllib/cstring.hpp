#ifndef BTLLIB_CSTRING_HPP
#define BTLLIB_CSTRING_HPP

#include <cstdlib>
#include <cstring>
#include <string>

namespace btllib {

struct CString
{
  static const size_t CSTRING_DEFAULT_CAP = 4096;

  CString() { s[0] = '\0'; }

  CString(const CString& cstring)
  {
    if (cstring.s_size > s_cap) {
      s = (char*)std::realloc((char*)s, cstring.s_size); // NOLINT
      s_cap = cstring.s_size;
    }
    s_size = cstring.s_size;
    memcpy(s, cstring.s, s_size);
  }

  CString(CString&& cstring) noexcept
  {
    std::swap(s, cstring.s);
    s_size = cstring.s_size;
    cstring.clear();
    std::swap(s_cap, cstring.s_cap);
  }

  CString(const std::string& str)
  {
    if (str.size() + 1 > s_cap) {
      s_cap = str.size() + 1;
      s = (char*)std::realloc((char*)s, s_cap); // NOLINT
    }
    s_size = str.size();
    memcpy(s, str.c_str(), s_size + 1);
  }

  CString& operator=(const CString& cstring)
  {
    if (this == &cstring) {
      return *this;
    }
    if (cstring.s_size > s_cap) {
      s = (char*)std::realloc((char*)s, cstring.s_size); // NOLINT
      s_cap = cstring.s_size;
    }
    s_size = cstring.s_size;
    memcpy(s, cstring.s, s_size);
    return *this;
  }

  CString& operator=(CString&& cstring) noexcept
  {
    std::swap(s, cstring.s);
    s_size = cstring.s_size;
    cstring.clear();
    std::swap(s_cap, cstring.s_cap);
    return *this;
  }

  CString& operator=(const std::string& str)
  {
    if (str.size() + 1 > s_cap) {
      s_cap = str.size() + 1;
      s = (char*)std::realloc((char*)s, s_cap); // NOLINT
    }
    s_size = str.size();
    memcpy(s, str.c_str(), s_size + 1);
    return *this;
  }

  ~CString() { free(s); } // NOLINT

  void clear()
  {
    s[0] = '\0';
    s_size = 0;
  }
  bool empty() const { return (ssize_t)s_size <= 0; }
  size_t size() const { return s_size; }

  operator char*() const { return s; }

  char* s = (char*)std::malloc(CSTRING_DEFAULT_CAP); // NOLINT
  size_t s_size = 0;
  size_t s_cap = CSTRING_DEFAULT_CAP;
};

} // namespace btllib

#endif