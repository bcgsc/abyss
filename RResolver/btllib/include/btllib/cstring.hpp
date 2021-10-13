#ifndef BTLLIB_CSTRING_HPP
#define BTLLIB_CSTRING_HPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
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
    std::memcpy(s, cstr.s, s_size + 1);
  }

  CString(CString&& cstr) noexcept
    : s_size(cstr.s_size)
  {
    std::swap(s, cstr.s);
    cstr.clear();
    std::swap(s_cap, cstr.s_cap);
  }

  explicit CString(const std::string& str)
  {
    if (str.size() + 1 > s_cap) {
      change_cap(str.size() + 1);
    }
    s_size = str.size();
    std::memcpy(s, str.c_str(), s_size + 1);
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
    std::memcpy(s, cstr.s, s_size + 1);
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
    std::memcpy(s, str.c_str(), s_size + 1);
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
    std::memcpy(s + s_size, cstr.s, cstr.s_size);
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
    std::memcpy(s + s_size, str.c_str(), str.size());
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
    s_size = 0;
    s[0] = '\0';
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

  void resize(const size_t n, const char c = '\0')
  {
    if (n > s_size) {
      change_cap(n + 1);
      for (size_t i = s_size; i < n; i++) {
        s[i] = c;
      }
    }
    s_size = n;
    s[s_size] = '\0';
  }

  CString& erase(const size_t pos = 0,
                 size_t len = std::numeric_limits<size_t>::max())
  {
    if (pos + len > size()) {
      len = size() - pos;
    }
    const ssize_t to_move = ssize_t(size()) - ssize_t(pos) - ssize_t(len);
    if (to_move > 0 && to_move < ssize_t(size())) {
      std::memmove(s + pos, s + pos + len, to_move);
    }
    resize(size() - len);
    return *this;
  }

  char* s = (char*)std::malloc(CSTRING_DEFAULT_CAP); // NOLINT
  size_t s_size = 0;
  size_t s_cap = CSTRING_DEFAULT_CAP;
};

} // namespace btllib

#endif