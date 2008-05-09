#ifndef PREPROCESSOPTIONS_H
#define PREPROCESSOPTIONS_H

#include <string>

namespace pp_opt {
	extern std::string outFile;
	extern std::string fastaFile;

	void parse(int argc, char* const* argv);
}

#endif
