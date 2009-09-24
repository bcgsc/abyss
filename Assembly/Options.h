#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstdio>
#include <string>
#include <vector>

namespace opt {
	extern int rank;
	extern int numProc;
	extern int kmerSize;
	extern unsigned erode;
	extern unsigned erodeStrand;
	extern int trimLen;
	extern float coverage;
	extern int bubbles;
	extern int chastityFilter;
	extern int trimMasked;
	extern std::string contigsPath;
	extern std::string contigsTempPath;
	extern std::string graphPath;
	extern std::string snpPath;
	extern FILE* snpFile;
	extern int verbose;
	extern std::vector<std::string> inFiles;
	extern bool colourSpace;

	void parse(int argc, char* const* argv);
}

#endif
