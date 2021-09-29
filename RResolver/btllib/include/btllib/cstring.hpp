#ifndef BTLLIB_CSTRING_HPP
#define BTLLIB_CSTRING_HPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

namespace btllib {

struct CString
{
  static const size_t CSTRING_DEFAULT_CAP = 2048;

  CString() { s[0] = '\0'; }

  CString(const CString& cstr)
  {
    if (cstr.s_size + 1 > s_cap) {
      change_cap(cstr.s_size + 1);
    }
    s_size = cstr.s_size;
    memcpy(s, cstr.s, s_size + 1);
  }

  CString(CString&& cstr) noexcept
  {
    std::swap(s, cstr.s);
    s_size = cstr.s_size;
    cstr.clear();
    std::swap(s_cap, cstr.s_cap);
  }

  CString(const std::string& str)
  {
    if (str.size() + 1 > s_cap) {
      change_cap(str.size() + 1);
    }
    s_size = str.size();
    memcpy(s, str.c_str(), s_size + 1);
  }

  CString& operator=(const CString& cstr)
  {
    if (this == &cstr) {
      return *this;
    }
    if (cstr.s_size + 1 > s_cap) {
      change_cap(cstr.s_size + 1);
    }
    s_size = cstr.s_size;
    memcpy(s, cstr.s, s_size + 1);
    return *this;
  }

  CString& operator=(CString&& cstr) noexcept
  {
    std::swap(s, cstr.s);
    s_size = cstr.s_size;
    cstr.clear();
    std::swap(s_cap, cstr.s_cap);
    return *this;
  }

  CString& operator=(const std::string& str)
  {
    if (str.size() + 1 > s_cap) {
      change_cap(str.size() + 1);
    }
    s_size = str.size();
    memcpy(s, str.c_str(), s_size + 1);
    return *this;
  }

  CString& operator+=(const CString& cstr)
  {
    const auto new_size = s_size + cstr.s_size;
    if (new_size + 1 > s_cap) {
      const auto factor =
        size_t(std::pow(2,
                        std::ceil(std::log2(double(new_size + 1)) -
                                  std::log2(double(s_size)))));
      change_cap(s_size * factor);
    }
    memcpy(s + s_size, cstr.s, cstr.s_size);
    s_size = new_size;
    return *this;
  }

  CString& operator+=(const std::string& str)
  {
    const auto new_size = s_size + str.size();
    if (new_size + 1 > s_cap) {
      const auto factor =
        size_t(std::pow(2,
                        std::ceil(std::log2(double(new_size + 1)) -
                                  std::log2(double(s_size)))));
      change_cap(s_size * factor);
    }
    memcpy(s + s_size, str.c_str(), str.size());
    s[new_size] = '\0';
    s_size = new_size;
    return *this;
  }

  CString& operator+=(const char c)
  {
    const auto new_size = s_size + 1;
    if (new_size + 1 > s_cap) {
      const auto factor =
        size_t(std::pow(2,
                        std::ceil(std::log2(double(new_size + 1)) -
                                  std::log2(double(s_size)))));
      change_cap(s_size * factor);
    }
    s[new_size - 1] = c;
    s[new_size] = '\0';
    s_size = new_size;
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
  char& front() const { return s[0]; }
  char& back() const { return s[s_size - 1]; }
  void pop_back()
  {
    s_size -= 1;
    s[s_size] = '\0';
  }

  operator char*() const { return s; }

  void change_cap(const size_t new_cap)
  {
    s_cap = new_cap;
    s = (char*)std::realloc(s, new_cap); // NOLINT
  }

  char* s = (char*)std::malloc(CSTRING_DEFAULT_CAP); // NOLINT
  size_t s_size = 0;
  size_t s_cap = CSTRING_DEFAULT_CAP;
};

} // namespace btllib

#endif