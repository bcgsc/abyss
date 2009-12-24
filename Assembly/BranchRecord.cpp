#include "BranchRecord.h"
#include <utility>

using namespace std;

/** Add a single sequence to the branch. */
void BranchRecord::addSequence(const PackedSeq& seq)
{
	m_data.push_back(seq);

	// Detect a loop by checking that the sequence is not already in the branch
	bool unique = m_seqMap.insert(seq.first).second;
	if(!unique)
		m_loopDetected = true;
}

/**
 * Remove all the sequences including and following the specified
 * iterator.
 */
void BranchRecord::truncate(BranchRecord::iterator position)
{
	ssize_t size = position - m_data.begin();
	assert(size >= 0);
	assert((size_t)size < m_data.size());
	(void)size;
	for (BranchData::const_iterator it = m_data.begin();
			it != m_data.end(); ++it)
		m_seqMap.erase(it->first);
	m_data.erase(position);
}

//
// Terminate the branch and indicate why
//
void BranchRecord::terminate(BranchState reason)
{
	assert(reason != BS_ACTIVE);
	m_state = reason;
}

/** Set the extensions and multiplicity of a sequence. */
void BranchRecord::setData(const PackedSeq& seqData)
{
	assert(m_data.back().first == seqData.first);
	m_data.back().second = seqData.second;
}

/** Forget the multiplicity information. */
void BranchRecord::clearMultiplicity()
{
	m_seqMap.clear();
}

//
// Get the branch length
//
size_t BranchRecord::getLength() const
{
	return m_data.size();
}

//
// Is the branch active?
//
bool BranchRecord::isActive() const
{
	return m_state == BS_ACTIVE;
}

//
// Get the state of the branch
//
BranchState BranchRecord::getState() const
{
	return m_state;	
}

//
// Get the direction of extension
//
extDirection BranchRecord::getDirection() const
{
	return m_dir;	
}

const Kmer& BranchRecord::getFirstSeq() const
{
	assert(!m_data.empty());
	return m_data.front().first;
}

//
// Get the last sequence in the data structure
// 
const Kmer& BranchRecord::getLastSeq() const
{
	assert(!m_data.empty());
	return m_data.back().first;
}

//
// Check if a sequence exists in the branch record
//
bool BranchRecord::exists(const Kmer& seq) const
{
	assert(!m_seqMap.empty());
	return m_seqMap.find(seq) != m_seqMap.end();
}

/** Check if the branch is too long. */
bool BranchRecord::isTooLong() const
{
	// If the maxLength == -1, no length check should be performed.
	return m_maxLength > -1 && getLength() > getMaxLength();
}

// Calculate the total multiplicity of this branch. The result is
// saved so that it may be fetched using getBranchMultiplicity after
// the multiplicity of each individual sequence has been forgotten.
// ignoreLast - this flag will be set when determining the branch 
// multiplicity for bubble removal. Since the multiplicity record
// lags behind the extension in a branch group, the last sequence
// will not have its multiplicity set when the branches join back together
// since that sequence is in all branches by definition it will not
// contribute any information to choosing a particular branch
// so it can be safely ignored
//
int BranchRecord::calculateBranchMultiplicity(bool ignoreLast)
{
	assert(!m_data.empty());
	BranchData::const_iterator endSeq = m_data.end();
	if(ignoreLast)
		endSeq--;

	int total = 0;
	for(BranchData::const_iterator iter = m_data.begin(); iter != endSeq; ++iter)
	{
		int m = iter->second.getMultiplicity();
		assert(m > 0);
		total += m;
	}
	assert(total > 0);
	return m_multiplicity = total;
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
	assert(getState() == BS_NOEXT || getState() == BS_LOOP);
	assert(getLength() > 1);
	Kmer first = getFirstSeq();
	Kmer last = getLastSeq();
	if (getDirection() == SENSE)
		last.reverseComplement();
	else
		first.reverseComplement();
	return first.compare(last) <= 0;
}
