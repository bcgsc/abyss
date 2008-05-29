#ifndef LOG_H
#define LOG_H

#include <fstream>

class Log
{
	public:
		Log(std::string filename);
		~Log();
		void write(std::string str);
		static void setID(int id) { m_id = id; }
		static int m_id;

	private:
		std::ofstream m_fileHandle;
};

int PrintDebug(int level, const char* format, ...);

#endif
