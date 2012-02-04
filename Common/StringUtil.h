#ifndef STRINGUTIL_H
#define STRINGUTIL_H 1

#include <cassert>
#include <iomanip>
#include <sstream>
#include <string>

/** Return the last character of s and remove it. */
static inline char chop(std::string& s)
{
	assert(s.length() > 1);
	unsigned back = s.length() - 1;
	char c = s[back];
	s.erase(back);
	return c;
}

/** If the last character of s is c, remove it and return true. */
static inline bool chomp(std::string& s, char c = '\n')
{
	unsigned back = s.length() - 1;
	if (!s.empty() && s[back] == c) {
		s.erase(back);
		return true;
	} else
		return false;
}

/** Return the SI representation of n. */
static inline std::string toSI(double n)
{
	std::ostringstream s;
	s << std::setprecision(3);
	if (n < 1e3)
		s << n << ' ';
	else if (n < 1e6)
		s << n/1e3 << " k";
	else if (n < 1e9)
		s << n/1e6 << " M";
	else if (n < 1e12)
		s << n/1e9 << " G";
	else
		s << n/1e12 << " T";
	return s.str();
}

/** Return the engineering string representation of n. */
template <typename T>
static inline std::string toEng(T n)
{
	std::ostringstream s;
	s << std::setprecision(4);
	if (n < 10000000)
		s << n;
	else if (n < 1e9)
		s << n/1e6 << "e6";
	else if (n < 1e12)
		s << n/1e9 << "e9";
	else
		s << n/1e12 << "e12";
	return s.str();
}

/** Return true if the second string is a prefix of the string s. */
template <size_t N>
bool startsWith(const std::string& s, const char (&prefix)[N])
{
	size_t n = N - 1;
	return s.size() > n && equal(s.begin(), s.begin() + n, prefix);
}

/** Return true if the second string is a suffix of the string s. */
template <size_t N>
bool endsWith(const std::string& s, const char (&suffix)[N])
{
	size_t n = N - 1;
	return s.size() > n && equal(s.end() - n, s.end(), suffix);
}

/** Return true if the second string is a suffix of the string s. */
static inline
bool endsWith(const std::string& s, const std::string& suffix)
{
	size_t n = suffix.size();
	return s.size() > n && equal(s.end() - n, s.end(),
			suffix.begin());
}

#endif
