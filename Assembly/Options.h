#ifndef ASSEMBLY_OPTIONS_H
#define ASSEMBLY_OPTIONS_H 1

#include <cstdio>
#include <string>
#include <vector>

namespace opt {
	extern int kmerSize;
	extern unsigned erode;
	extern unsigned erodeStrand;
	extern int trimLen;
	extern float coverage;
	extern int bubbles;
	extern std::string coverageHistPath;
	extern std::string contigsPath;
	extern std::string contigsTempPath;
	extern std::string graphPath;
	extern std::string snpPath;
	extern FILE* snpFile;
	extern std::vector<std::string> inFiles;

	void parse(int argc, char* const* argv);
}

#endif
