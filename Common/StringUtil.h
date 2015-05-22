#ifndef STRINGUTIL_H
#define STRINGUTIL_H 1

#include <cassert>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

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

/** Return the SI representation of a number in bytes. */
static inline std::string bytesToSI(size_t n)
{
	std::ostringstream s;
	s << std::setprecision(3);
	if (n < 1024)
		s << n;
	else if (n < (1ULL<<20))
		s << (double)n/(1ULL<<10) << "k";
	else if (n < (1ULL<<30))
		s << (double)n/(1ULL<<20) << "M";
	else
		s << (double)n/(1ULL<<30) << "G";
	return s.str();
}

/**
 * Convert a quantity with SI units to the equivalent floating
 * point number.
 */
static inline double fromSI(std::istringstream& iss)
{
	double size;
	std::string units;

	iss >> size;
	if (iss.fail()) {
		// not prefixed by a number
		return 0;
	}

	iss >> units;
	if (iss.fail() && iss.eof()) {
		// no units given; clear fail flag
		// and just return the number
		iss.clear(std::ios::eofbit);
		return size;
	}

	if (units.size() > 1) {
		// unrecognized multichar suffix
		iss.setstate(std::ios::failbit);
		return 0;
	}

	switch(tolower(units[0])) {
		case 'k':
			size *= 1000ULL; break;
		case 'm':
			size *= 1000ULL*1000; break;
		case 'g':
			size *= 1000ULL*1000*1000; break;
		case 't':
			size *= 1000ULL*1000*1000*1000; break;
		default:
			iss.setstate(std::ios::failbit);
			return 0;
	}

	return size;
}

/**
 * Convert a quantity with SI units to the equivalent floating
 * point number.
 */
static inline double fromSI(const std::string& str)
{
	std::istringstream iss(str);
	return fromSI(iss);
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

static inline
bool isReadNamePair(const std::string& name1, const std::string& name2)
{
	assert(!name1.empty() && !name2.empty());

	if (name1 == name2)
		return true;

	if (endsWith(name1,"/1") && endsWith(name2,"/2")) {
		int len1 = name1.length();
		int len2 = name2.length();
		assert(len1 > 2 && len2 > 2);
		return name1.compare(0, len1-2, name2, 0, len2-2) == 0;
	}

	return false;
}

static inline size_t SIToBytes(std::istringstream& iss)
{
	double size;
	std::string units;

	iss >> size;
	if (iss.fail()) {
		// not prefixed by a number
		return 0;
	}

	iss >> units;
	if (iss.fail() && iss.eof()) {
		// no units given; clear fail flag
		// and assume bytes
		iss.clear(std::ios::eofbit);
		return (size_t)ceil(size);
	}

	if (units.size() > 1) {
		// unrecognized multichar suffix
		iss.setstate(std::ios::failbit);
		return 0;
	}

	switch(tolower(units[0])) {
		case 'k':
			size *= (size_t)1<<10; break;
		case 'm':
			size *= (size_t)1<<20; break;
		case 'g':
			size *= (size_t)1<<30; break;
		default:
			iss.setstate(std::ios::failbit);
			return 0;
	}

	return (size_t)ceil(size);
}

static inline size_t SIToBytes(const std::string& str)
{
	std::istringstream iss(str);
	return SIToBytes(iss);
}

#endif
