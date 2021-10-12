/**
 * Random utility functions.
 */
#ifndef BTLLIB_UTIL_HPP
#define BTLLIB_UTIL_HPP

#include "cstring.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace btllib {

inline std::vector<std::string>
split(const std::string& s, const std::string& delim);
inline std::string
join(const std::vector<std::string>& s, const std::string& delim);
inline void
ltrim(std::string& s);
inline void
ltrim(CString& s);
inline void
rtrim(std::string& s);
inline void
rtrim(CString& s);
inline void
trim(std::string& s);
inline void
trim(CString& s);
inline bool
startswith(std::string s, std::string prefix);
inline bool
endswith(std::string s, std::string suffix);

inline std::vector<std::string>
split(const std::string& s, const std::string& delim)
{
  std::vector<std::string> tokens;
  size_t pos1 = 0, pos2 = 0;
  while ((pos2 = s.find(delim, pos2)) != std::string::npos) {
    tokens.push_back(s.substr(pos1, pos2 - pos1));
    pos2 += delim.size();
    pos1 = pos2;
  }
  tokens.push_back(s.substr(pos1));
  return tokens;
}

inline std::string
join(const std::vector<std::string>& s, const std::string& delim)
{
  std::string joined = s[0];
  for (size_t i = 1; i < s.size(); i++) {
    joined += delim;
    joined += s[i];
  }
  return joined;
}

inline void
ltrim(std::string& s)
{
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
            return !bool(std::isspace(ch));
          }));
}

inline void
ltrim(CString& s)
{
  decltype(s.size()) i = 0;
  while (i < s.size() && bool(std::isspace(s[i]))) {
    i++;
  }
  s.erase(0, i);
}

inline void
rtrim(std::string& s)
{
  s.erase(std::find_if(s.rbegin(),
                       s.rend(),
                       [](int ch) { return !bool(std::isspace(ch)); })
            .base(),
          s.end());
}

inline void
rtrim(CString& s)
{
  auto i = s.size();
  while (i > 0 && bool(std::isspace(s[i - 1]))) {
    i--;
  }
  s.resize(i);
}

inline void
trim(std::string& s)
{
  ltrim(s);
  rtrim(s);
}

inline void
trim(CString& s)
{
  ltrim(s);
  rtrim(s);
}

inline bool
startswith(std::string s, std::string prefix)
{
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  std::transform(prefix.begin(), prefix.end(), prefix.begin(), ::tolower);
  return s.find(prefix) == 0;
}

inline bool
endswith(std::string s, std::string suffix)
{
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  std::transform(suffix.begin(), suffix.end(), suffix.begin(), ::tolower);
  auto pos = s.rfind(suffix);
  return (pos != std::string::npos) && (pos == s.size() - suffix.size());
}

} // namespace btllib

#endif