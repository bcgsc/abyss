#ifndef BTLLIB_CSTRING_HPP
#define BTLLIB_CSTRING_HPP

#include <cstdlib>
#include <limits>
#include <string>

namespace btllib {

struct CString
{
  static const size_t CSTRING_DEFAULT_CAP = 2048;

  CString() { s[0] = '\0'; }
  CString(const CString& cstr);
  CString(CString&& cstr) noexcept;
  explicit CString(const std::string& str);

  CString& operator=(const CString& cstr);
  CString& operator=(CString&& cstr) noexcept;
  CString& operator=(const std::string& str);
  CString& operator+=(const CString& cstr);
  CString& operator+=(const std::string& str);
  CString& operator+=(char c);

  ~CString() { free(s); } // NOLINT

  void clear();
  bool empty() const { return (ssize_t)s_size <= 0; }
  size_t size() const { return s_size; }
  char& front() const { return s[0]; }
  char& back() const { return s[s_size - 1]; }
  void pop_back();

  operator char*() const { return s; }

  void change_cap(size_t new_cap);
  void resize(size_t n, char c = '\0');

  CString& erase(size_t pos = 0,
                 size_t len = std::numeric_limits<size_t>::max());

  char* s = (char*)std::malloc(CSTRING_DEFAULT_CAP); // NOLINT
  size_t s_size = 0;
  size_t s_cap = CSTRING_DEFAULT_CAP;
};

} // namespace btllib

#endif