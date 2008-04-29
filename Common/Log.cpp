#include "Log.h"

Log::Log(std::string filename)
{
	std::ios_base::openmode mode = std::ios::out;
	m_fileHandle.open(filename.c_str(), mode);
	assert(m_fileHandle.is_open());	
}

Log::~Log()
{
	m_fileHandle.close();
}

void Log::write(std::string str)
{
	m_fileHandle << str << std::endl;
	
	// write out the message immediately
	m_fileHandle.flush();
}
