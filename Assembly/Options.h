#ifndef ASSEMBLY_OPTIONS_H
#define ASSEMBLY_OPTIONS_H 1

#include <string>
#include <vector>

namespace opt {
	extern unsigned kmerSize;
	extern unsigned singleKmerSize;
	extern unsigned kMin;
	extern unsigned kMax;
	extern unsigned kStep;
	extern unsigned erode;
	extern unsigned erodeStrand;
	extern unsigned trimLen;
	extern float coverage;
	extern unsigned kc;
	extern unsigned bubbleLen;
	extern unsigned ss;
	extern bool maskCov;
	extern std::string coverageHistPath;
	extern std::string contigsPath;
	extern std::string contigsTempPath;
	extern std::string graphPath;
	extern std::string snpPath;
	extern std::vector<std::string> inFiles;

	extern std::string db;

	void parse(int argc, char* const* argv);
	extern std::string assemblyCmd;
	std::vector<std::string> getMetaValue();
	int getVvalue();
	std::string getUvalue();
	std::string getCommand();
}

#endif
