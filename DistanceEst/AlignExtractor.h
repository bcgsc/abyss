#ifndef ALIGNEXTRACTOR
#define ALIGNEXTRACTOR 1

#include "SAM.h"
#include <istream>
#include <ostream>
#include <string>

typedef SAMRecord AlignPair;

typedef std::vector<AlignPair> AlignPairVec;

class AlignExtractor
{
	public:
		AlignExtractor(std::istream& in);

		// return true if EOF was reached
		bool extractContigAlignments(AlignPairVec& outPairs);

	private:
		std::istream& m_in;
		AlignPair m_currPair;
};

#endif
