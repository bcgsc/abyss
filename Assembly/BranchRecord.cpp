#include "BranchRecord.h"
#include <utility>

using namespace std;

/** Return true if any of the last three k-mer of this branch is the
 * specified kmer. */
bool BranchRecord::exists(const Kmer& kmer) const
{
	assert(!m_data.empty());
	return m_data.back().first == kmer
		|| (m_data.size() > 1 && (m_data.rbegin()+1)->first == kmer)
		|| (m_data.size() > 2 && (m_data.rbegin()+2)->first == kmer);
}

/** Check if the branch is too long. */
bool BranchRecord::isTooLong() const
{
	// If the maxLength == -1, no length check should be performed.
	return m_maxLength > -1 && size() > getMaxLength();
}

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
	return first.compare(last) <= 0;
}
