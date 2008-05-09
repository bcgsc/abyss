#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>

namespace opt {
	extern int kmerSize;
	extern int readLen;
	extern int trimLen;
	extern bool disableErosion;
	extern std::string inFile;

	void parse(int argc, char* const* argv);
}

#endif
