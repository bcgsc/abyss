#ifndef ALIGNEXTRACtOR
#define ALIGNEXTRACTOR

#include <iostream>
#include <fstream>
#include "AlignmentCache.h"

struct AlignPair
{
	AlignData refRec;
	AlignData pairRec;
};

typedef std::vector<AlignPair> AlignPairVec;


class AlignExtractor
{
	public:
		AlignExtractor(std::string file);
		~AlignExtractor();

		// return true if EOF was reached
		bool extractContigAlignments(AlignPairVec& outPairs);

	private:
		AlignPair readRecord();

		AlignPair m_currPair;
		std::ifstream m_fileHandle;

};

#endif
