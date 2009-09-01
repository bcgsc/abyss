#include "Timer.h"
#include "Log.h"

using namespace std;

// Constructor starts the timer
Timer::Timer(string funcString)
	: m_funcStr(funcString), m_start(clock())
{
}

// Destructor stops it and prints
Timer::~Timer()
{
	PrintDebug(2, "%s: %.3f s\n", m_funcStr.c_str(),
			(double)(clock() - m_start) / CLOCKS_PER_SEC);
}
