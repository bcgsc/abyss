#include "Log.h"
#include <stdarg.h>
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

int Log::m_id = -1;

int PrintDebug(int level, const char* format, ...)
{
	if (level > 3)
		return 0;
	if (Log::m_id >= 0)
		printf("%d: ", Log::m_id);
	va_list ap;
	va_start(ap, format);
	int retval = vprintf(format, ap);
	va_end(ap);
	return retval;
}
