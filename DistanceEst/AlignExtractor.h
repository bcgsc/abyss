#ifndef ALIGNEXTRACTOR
#define ALIGNEXTRACTOR 1

#include "SAM.h"
#include <istream>

typedef std::vector<SAMRecord> AlignPairVec;

class AlignExtractor
{
	public:
		AlignExtractor(std::istream& in);

		bool extractContigAlignments(AlignPairVec& outPairs);

		operator void*() { return m_in; }

		friend AlignExtractor& operator >>(
				AlignExtractor& in, AlignPairVec& outPairs)
		{
			in.extractContigAlignments(outPairs);
			return in;
		}

	private:
		std::istream& m_in;
		SAMRecord m_currPair;
};

#endif
