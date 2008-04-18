#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>

namespace opt {
	extern unsigned kmerSize;
	extern unsigned readLen;
	extern std::string fastaFile;
	void parse(int argc, char* const* argv);
}

#endif
