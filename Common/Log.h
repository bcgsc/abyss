#ifndef LOG_H
#define LOG_H

#include <fstream>

class Log
{
	public:
		Log(std::string filename);
		~Log();

		void write(std::string str);

	private:
		std::ofstream m_fileHandle;
};

#endif
