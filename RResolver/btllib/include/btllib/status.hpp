/**
 * Functions for logging and error checking.
 */
#ifndef BTLLIB_STATUS_HPP
#define BTLLIB_STATUS_HPP

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <string>

namespace btllib {

constexpr const char* PRINT_COLOR_INFO = "\33[32m";
constexpr const char* PRINT_COLOR_WARNING = "\33[33m";
constexpr const char* PRINT_COLOR_ERROR = "\33[31m";
constexpr const char* PRINT_COLOR_END = "\33[0m";

inline std::string
get_time();

/**
 * Log info level events.
 *
 * @param msg Message to print.
 */
inline void
log_info(const std::string& msg);

/**
 * Log warning level events.
 *
 * @param msg Message to print.
 */
inline void
log_warning(const std::string& msg);

/**
 * Log error level events.
 *
 * @param msg Message to print.
 */
inline void
log_error(const std::string& msg);

/**
 * Conditionally log info level events.
 *
 * @param condition If this is true, the message is printed.
 * @param msg Message to print.
 */
inline void
check_info(bool condition, const std::string& msg);

/**
 * Conditionally log warning level events.
 *
 * @param condition If this is true, the message is printed.
 * @param msg Message to print.
 */
inline void
check_warning(bool condition, const std::string& msg);

/**
 * Conditionally log error level events. The program exits if the condition is
 * true.
 *
 * @param condition If this is true, the message is printed and the program
 * exits.
 * @param msg Message to print.
 */
inline void
check_error(bool condition, const std::string& msg);

/**
 * Check whether the stream is good. Program prints an error message and exits
 * if not.
 *
 * @param stream Stream to check goodness of.
 * @param name Name of the stream, e.g. filepath or stdin
 */
inline void
check_stream(const std::ios& stream, const std::string& name);

inline std::string
get_time()
{
  time_t now;
  time(&now);
  char buf[sizeof("2011-10-08T07:07:09Z")];
  std::tm tm_result = {};
  localtime_r(&now, &tm_result);
  std::strftime(buf, sizeof buf, "%F %T", &tm_result);
  return std::string(buf);
}

inline void
log_info(const std::string& msg)
{
  std::cerr << ('[' + get_time() + "]" + PRINT_COLOR_INFO + "[INFO] " +
                PRINT_COLOR_END + msg + '\n')
            << std::flush;
}

inline void
log_warning(const std::string& msg)
{
  std::cerr << ('[' + get_time() + "]" + PRINT_COLOR_WARNING + "[WARNING] " +
                PRINT_COLOR_END + msg + '\n')
            << std::flush;
}

inline void
log_error(const std::string& msg)
{
  std::cerr << ('[' + get_time() + "]" + PRINT_COLOR_ERROR + "[ERROR] " +
                PRINT_COLOR_END + msg + '\n')
            << std::flush;
}

inline void
check_info(bool condition, const std::string& msg)
{
  if (condition) {
    log_info(msg);
  }
}

inline void
check_warning(bool condition, const std::string& msg)
{
  if (condition) {
    log_warning(msg);
  }
}

inline void
check_error(bool condition, const std::string& msg)
{
  if (condition) {
    log_error(msg);
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
}

inline std::string
get_strerror()
{
  static const size_t buflen = 1024;
  char buf[buflen];
// POSIX and GNU implementation of strerror_r differ, even in function signature
// and so we need to check which one is used
#if __APPLE__ ||                                                               \
  ((_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && !_GNU_SOURCE)
  strerror_r(errno, buf, buflen);
  return buf;
#else
  return strerror_r(errno, buf, buflen);
#endif
}

inline void
check_stream(const std::ios& stream, const std::string& name)
{
  if (!stream.good()) {
    log_error("'" + name + "' stream error: " + get_strerror());
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
}

} // namespace btllib

#endif
