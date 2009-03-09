#ifndef TIMER_H
#define TIMER_H

#include <string>

// Simple class to automatically time how long a function takes

class Timer
{
	public:
		Timer(std::string funcString);
		~Timer();
		
	private:
		std::string m_funcStr;
		clock_t m_start;
};

#endif
