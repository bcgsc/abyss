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

/** Build a contig from a branch. */
BranchRecord::operator Sequence() const
{
	assert(!m_data.empty());
	Sequence outseq;
	size_t n = m_data.front().first.length() + m_data.size() - 1;
	outseq.reserve(n);

	if (m_dir == SENSE) {
		BranchData::const_iterator iter = m_data.begin();
		outseq = iter->first.str();
		outseq.resize(n, 'N');
		++iter;
		string::iterator outa = outseq.begin() + Kmer::length();
		string::iterator outb = outseq.begin() + KmerPair::length();
		for (; iter != m_data.end(); ++iter) {
			std::pair<char, char> x = iter->first.getLastBaseChar();
			assert(*outa == 'N' || *outa == x.first);
			*outa = x.first;
			++outa;
			assert(outb < outseq.end());
			assert(*outb == 'N');
			*outb = x.second;
			++outb;
		}
	} else {
		BranchData::const_reverse_iterator iter = m_data.rbegin();
		outseq = iter->first.str();
		outseq.resize(n, 'N');
		++iter;
		string::iterator outa = outseq.begin() + Kmer::length();
		string::iterator outb = outseq.begin() + KmerPair::length();
		for (; iter != m_data.rend(); ++iter) {
			std::pair<char, char> x = iter->first.getLastBaseChar();
			assert(*outa == 'N' || *outa == x.first);
			*outa = x.first;
			++outa;
			assert(outb < outseq.end());
			assert(*outb == 'N');
			*outb = x.second;
			++outb;
		}
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
	V first = front().first;
	V last = back().first;
	if (getDirection() == SENSE)
		last.reverseComplement();
	else
		first.reverseComplement();
	assert(first != last);
	return first < last;
}
