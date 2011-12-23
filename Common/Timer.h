#ifndef TIMER_H
#define TIMER_H 1

#include <string>

/**
 * Time the duration between the construction and destruction of this
 * timer object and log that duration.
 */
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
