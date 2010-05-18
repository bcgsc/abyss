#include "AlignExtractor.h"
#include <cassert>

using namespace std;

AlignExtractor::AlignExtractor(istream& in)
	: m_in(in)
{
	// Ignore SAM headers.
	while (m_in.peek() == '@') {
		m_in.ignore(numeric_limits<streamsize>::max(), '\n');
		assert(m_in.good());
	}

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
	if (m_in.eof())
		m_in.setstate(ios_base::failbit);
	if (!m_in)
		return true;
	outPairs.clear();
	outPairs.push_back(m_currPair);
	ContigID id = m_currPair.rname;
	for (SAMRecord pair; m_in >> pair;) {
		if (pair.rname == id) {
			outPairs.push_back(pair);
		} else {
			m_currPair = pair;
			return false;
		}
	}
	assert(m_in.eof());
	m_in.clear(ios_base::eofbit);
	return true;
}
