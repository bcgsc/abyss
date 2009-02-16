#ifndef ALIGNEXTRACTOR
#define ALIGNEXTRACTOR 1

#include "Aligner.h"
#include <istream>
#include <ostream>
#include <string>

struct AlignPair
{
	Alignment refRec;
	Alignment pairRec;
	friend std::istream& operator >>(std::istream& in, AlignPair& p)
	{
		std::string ida, idb;
		return in >> ida >> p.refRec >> idb >> p.pairRec;
	}
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
		AlignExtractor(std::istream& in);

		// return true if EOF was reached
		bool extractContigAlignments(AlignPairVec& outPairs);

	private:
		std::istream& m_in;
		AlignPair m_currPair;
};

#endif
