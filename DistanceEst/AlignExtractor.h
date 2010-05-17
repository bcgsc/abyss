#ifndef ALIGNEXTRACTOR
#define ALIGNEXTRACTOR 1

#include "SAM.h"
#include <istream>

typedef std::vector<SAMRecord> AlignPairVec;

class AlignExtractor
{
	public:
		AlignExtractor(std::istream& in);

		// return true if EOF was reached
		bool extractContigAlignments(AlignPairVec& outPairs);

	private:
		std::istream& m_in;
		SAMRecord m_currPair;
};

#endif
