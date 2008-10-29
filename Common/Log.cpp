#include "config.h"
#include "Log.h"
#include "Options.h"
#include <cassert>
#include <cstdarg>
#include <unistd.h>

Log::Log(std::string filename)
	: m_fileHandle(filename.c_str())
{
	assert(m_fileHandle.is_open());
	char hostname[HOST_NAME_MAX];
	gethostname(hostname, sizeof hostname);
	m_fileHandle << hostname << std::endl;
}

Log::~Log()
{
	m_fileHandle.close();
}

void Log::write(std::string str)
{
	m_fileHandle << str << std::endl;
}

int PrintDebug(int level, const char* format, ...)
{
	if (opt::verbose < level)
		return 0;
	if (opt::rank >= 0)
		printf("%d: ", opt::rank);
	va_list ap;
	va_start(ap, format);
	int retval = vprintf(format, ap);
	va_end(ap);
	return retval;
}
