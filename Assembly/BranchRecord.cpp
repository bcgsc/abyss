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
		assert(m > 0);
		total += m;
	}
	assert(total > 0);
	return total;
}

/** Build a contig from a branch. */
BranchRecord::operator Sequence() const
{
	assert(!m_data.empty());
	Sequence outseq;
	outseq.reserve(m_data.front().first.length() + m_data.size() - 1);

	if (m_dir == SENSE) {
		BranchData::const_iterator iter = m_data.begin();
		outseq = iter->first.decode();
		++iter;
		for (; iter != m_data.end(); ++iter)
			outseq.append(1, iter->first.getLastBaseChar());
	} else {
		BranchData::const_reverse_iterator iter = m_data.rbegin();
		outseq = iter->first.decode();
		++iter;
		for (; iter != m_data.rend(); ++iter)
			outseq.append(1, iter->first.getLastBaseChar());
	}
	return outseq;
}

/**
 * Return whether this branch is the canonical representation of the
 * contig that it represents. A contig has two ends, and the contig
 * is output starting from the lexicographically smaller end.
 */
bool BranchRecord::isCanonical() const
{
	assert(size() > 1);
	Kmer first = front().first;
	Kmer last = back().first;
	if (getDirection() == SENSE)
		last.reverseComplement();
	else
		first.reverseComplement();
	assert(first != last);
	return first < last;
}
