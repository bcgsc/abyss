#ifndef ALIGNEXTRACTOR
#define ALIGNEXTRACTOR 1

#include <iostream>
#include <fstream>
#include "AlignmentCache.h"

struct AlignPair
{
	Alignment refRec;
	Alignment pairRec;
	friend std::ostream& operator <<(std::ostream& o,
			const AlignPair& p)
	{
		return o << p.refRec << ' ' << p.pairRec;
	}
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
