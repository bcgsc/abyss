/**
 * Random utility functions.
 */
#ifndef BTLLIB_UTIL_HPP
#define BTLLIB_UTIL_HPP

#include "btllib/cstring.hpp"

#include <condition_variable>
#include <mutex>
#include <string>
#include <vector>

namespace btllib {

/**
 * Split a string into component substrings with a delimiter.
 *
 * @param s String to split.
 * @param delim Delimiter to split with.
 *
 * @return Vector of substrings delimited by `delim`, excluding delimiters
 * themselves.
 */
std::vector<std::string>
split(const std::string& s, const std::string& delim);

/**
 * Join a vector of strings into a single string with a delimiter.
 *
 * @param s Vector of strings to join.
 * @param delim Delimiter to join the strings with.
 *
 * @return String with all the components joined.
 */
std::string
join(const std::vector<std::string>& s, const std::string& delim);

/**
 * Trim whitespace on the left side of the given string.
 *
 * @param s String to trim, edited in-place.
 *
 */
void
ltrim(std::string& s);
void
ltrim(btllib::CString& s);

/**
 * Trim whitespace on the right side of the given string.
 *
 * @param s String to trim, edited in-place.
 *
 */
void
rtrim(std::string& s);
void
rtrim(btllib::CString& s);

/**
 * Trim whitespace on the left and right side of the given string.
 *
 * @param s String to trim, edited in-place.
 *
 */
void
trim(std::string& s);
void
trim(btllib::CString& s);

/**
 * Check whether the given string starts with a prefix.
 *
 * @param s String to check.
 * @param prefix Prefix to check for.
 *
 */
bool
startswith(std::string s, std::string prefix);

/**
 * Check whether the given string ends with a suffix.
 *
 * @param s String to check.
 * @param suffix Suffix to check for.
 *
 */
bool
endswith(std::string s, std::string suffix);

/**
 * Equivalent to the GNU implementation of basename,
 * but returns a string copy of the result.
 *
 * @param path The path to get basename from.
 *
 * @return The basename of the path.
 */
std::string
get_basename(const std::string& path);

/**
 * Equivalent to the GNU implementation of dirname,
 * but returns a string copy of the result.
 *
 * @param path The path to get dirname from.
 *
 * @return The dirname of the path.
 */
std::string
get_dirname(const std::string& path);

// This exists in C++20, but we don't support that yet
/// @cond HIDDEN_SYMBOLS
class Barrier
{

public:
  Barrier(const unsigned count)
    : counter_default(count)
  {
  }

  void wait();

private:
  std::mutex m;
  std::condition_variable cv;
  unsigned counter{ 0 };
  unsigned counter_default;
  unsigned waiting{ 0 };
};
/// @endcond

} // namespace btllib

#endif