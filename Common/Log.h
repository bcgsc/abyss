#ifndef LOG_H
#define LOG_H 1

int PrintDebug(int level, const char* format, ...)
	__attribute__((format(printf, 2, 3)));

#endif
