#ifndef LOG_H
#define LOG_H 1

#include <ostream>

std::ostream& clog(int level);

int PrintDebug(int level, const char* format, ...)
	__attribute__((format(printf, 2, 3)));

#endif
