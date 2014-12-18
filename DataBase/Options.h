#ifndef OPTIONS_H
#define OPTIONS_H 1

#include <string>
#include <sstream>
#include <iterator>

namespace opt {
	extern std::string db;

	std::string getCommand(
			int argc, char* const* argv)
	{
		std::ostringstream command;
		char* const* last = argv + argc -1;
		copy(argv, last, std::ostream_iterator<const char*>(command, " "));
		command << *last;
		return command.str();
	}
}

#endif
