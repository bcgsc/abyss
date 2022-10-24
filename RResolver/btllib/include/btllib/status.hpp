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

/**
 * Return current time as a string.
 */
std::string
get_time();

/**
 * Log info level events.
 *
 * @param msg Message to print.
 */
void
log_info(const std::string& msg);

/**
 * Log warning level events.
 *
 * @param msg Message to print.
 */
void
log_warning(const std::string& msg);

/**
 * Log error level events.
 *
 * @param msg Message to print.
 */
void
log_error(const std::string& msg);

/**
 * Conditionally log info level events.
 *
 * @param condition If this is true, the message is printed.
 * @param msg Message to print.
 */
void
check_info(bool condition, const std::string& msg);

/**
 * Conditionally log warning level events.
 *
 * @param condition If this is true, the message is printed.
 * @param msg Message to print.
 */
void
check_warning(bool condition, const std::string& msg);

/**
 * Conditionally log error level events. The program exits if the condition is
 * true.
 *
 * @param condition If this is true, the message is printed and the program
 * exits.
 * @param msg Message to print.
 */
void
check_error(bool condition, const std::string& msg);

std::string
get_strerror();

/**
 * Check whether the stream is good. Program prints an error message and exits
 * if not.
 *
 * @param stream Stream to check goodness of.
 * @param name Name of the stream, e.g. filepath or stdin
 */
void
check_stream(const std::ios& stream, const std::string& name);

/**
 * Check whether the file at the given path is accessible (exists, permissions
 * are good, etc.).
 *
 * @param filepath Path to the file to check.
 */
void
check_file_accessibility(const std::string& filepath);

} // namespace btllib

#endif
