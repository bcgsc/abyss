#include "Log.h"
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
