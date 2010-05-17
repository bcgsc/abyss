#include "AlignExtractor.h"
#include <cassert>

using namespace std;

AlignExtractor::AlignExtractor(istream& in)
	: m_in(in)
{
	// Prime the read by reading in the first contig.
	bool good = m_in >> m_currPair;
	assert(good);
	(void)good;
}

/** Read alignment pairs and store them in the specified vector.
 * @return true at end-of-file
 */
bool AlignExtractor::extractContigAlignments(AlignPairVec& outPairs)
{
	assert(m_in.good());
	outPairs.push_back(m_currPair);
	ContigID id = m_currPair.rname;
	AlignPair pair;
	while (m_in >> pair) {
		if (pair.rname == id) {
			outPairs.push_back(pair);
		} else {
			m_currPair = pair;
			return false;
		}
	}
	assert(m_in.eof());
	return true;
}
