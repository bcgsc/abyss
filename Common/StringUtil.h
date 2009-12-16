#ifndef STRINGUTIL_H
#define STRINGUTIL_H 1

#include <cassert>
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

#endif
