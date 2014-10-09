#include "BranchRecord.h"

using namespace std;

/** Calculate the total multiplicity of this branch. */
int BranchRecord::calculateBranchMultiplicity() const
{
	assert(!m_data.empty());
	int total = 0;
	for (BranchData::const_iterator it = m_data.begin();
			it != m_data.end(); ++it) {
		int m = it->second.getMultiplicity();
		assert(m >= 0);
		total += m;
	}
	assert(total >= 0);
	return total;
}

/**
 * Return whether this branch is the canonical representation of the
 * contig that it represents. A contig has two ends, and the contig
 * is output starting from the lexicographically smaller end.
 */
bool BranchRecord::isCanonical() const
{
	assert(size() > 1);
	V first = front().first;
	V last = back().first;
	if (getDirection() == SENSE)
		last.reverseComplement();
	else
		first.reverseComplement();
	assert(first != last);
	return first < last;
}
