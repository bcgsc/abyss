#include "BranchRecord.h"
#include <iostream>
//
// Constructor
//
BranchRecord::BranchRecord()
	: m_dir(SENSE), m_state(BS_ACTIVE), m_maxLength(-1),
	m_loopDetected(false), m_multiplicity(-1)
{
}

//
// Constructor
//
BranchRecord::BranchRecord(extDirection dir, int maxLength)
	: m_dir(dir), m_state(BS_ACTIVE), m_maxLength(maxLength),
	m_loopDetected(false), m_multiplicity(-1)
{
}

//
// Copy Constructor
//
BranchRecord::BranchRecord(const BranchRecord& other)
{
	m_data = other.m_data;
	m_seqMap = other.m_seqMap;
	m_dir = other.m_dir;
	m_maxLength = other.m_maxLength;
	m_state = other.m_state;
	m_loopDetected = other.m_loopDetected;
	m_multiplicity = other.m_multiplicity;
}

// Assignment operator
BranchRecord& BranchRecord::operator=(const BranchRecord& other)
{
	// Detect self assignment
	if (this == &other)
	{
		return *this;
	}
		
	m_data = other.m_data;
	m_seqMap = other.m_seqMap;
	m_dir = other.m_dir;
	m_maxLength = other.m_maxLength;
	m_state = other.m_state;
	m_loopDetected = other.m_loopDetected;	
	m_multiplicity = other.m_multiplicity;
	return *this;
}

//
// Add a single sequence to the branch
//
void BranchRecord::addSequence(const PackedSeq& seq, int multiplicity)
{
	m_data.push_back(seq);
	
	// Detect a loop by checking that the sequence is not already in the branch
	MultMapPair item(seq, multiplicity);
	bool unique = m_seqMap.insert(item).second;
	if(!unique)
	{
		m_loopDetected = true;
	}
}

/**
 * Remove all the sequences including and following the specified
 * iterator.
 */
void BranchRecord::truncate(BranchDataIter position)
{
	ssize_t size = position - m_data.begin();
	assert(size > 0);
	assert((size_t)size < m_data.size());
	m_data.resize(size);
}

//
// Terminate the branch and indicate why
//
void BranchRecord::terminate(BranchState reason)
{
	assert(reason != BS_ACTIVE);
	m_state = reason;
}

//
// Get the multiplicity of the sequence
//
int BranchRecord::getMultiplicity(const PackedSeq& seq) const
{
	BranchMultMap::const_iterator iter = m_seqMap.find(seq);
	assert(iter != m_seqMap.end());
	return iter->second;
}

//
// Set the multiplicity of a sequence
//
void BranchRecord::setMultiplicity(const PackedSeq& seq, int multiplicity)
{
	BranchMultMap::iterator iter = m_seqMap.find(seq);
	assert(iter != m_seqMap.end());
	iter->second = multiplicity;
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

const PackedSeq& BranchRecord::getFirstSeq() const
{
	assert(!m_data.empty());
	return m_data.front();	
}

//
// Get the last sequence in the data structure
// 
const PackedSeq& BranchRecord::getLastSeq() const
{
	assert(!m_data.empty());
	return m_data.back();
}

//
// Check if a sequence exists in the branch record
//
bool BranchRecord::exists(const PackedSeq& seq) const
{
	assert(!m_seqMap.empty());
	return m_seqMap.find(seq) != m_seqMap.end();
}

//
// If the maxLength == -1, no length check should be performed
//
bool BranchRecord::doLengthCheck() const
{
	return (m_maxLength > -1);
}

/** Check if the branch is too long. */
bool BranchRecord::isTooLong() const
{
	return doLengthCheck() && getLength() > getMaxLength();
}

//
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
	assert(!m_seqMap.empty());
	int total = 0;
	
	BranchData::const_iterator endSeq = m_data.end();
	
	if(ignoreLast)
	{
		endSeq--;
	}
	
	size_t idx = 0;
	size_t numKmers = m_data.size();
	
	for(BranchData::const_iterator iter = m_data.begin(); iter != endSeq; ++iter)
	{
		int m = getMultiplicity(*iter);
		if(m <= 0)
		{
			std::cerr << "Multipliticty check failed! Node idx: " << idx << " branch size(kmers): " << numKmers << "\n";
			std::cerr << "Ignore last flag: " << ignoreLast << " BranchState: " << m_state << "\n";
			std::cerr << "Branch: ";
			printBranch(std::cerr);
			std::cerr << "\n";
			assert(m > 0);
		}
		
		total += m;
		++idx;
	}
	
	assert(total > 0);
	return m_multiplicity = total;
}

/** Build a contig from a branch. */
void BranchRecord::buildContig(Sequence& outseq) const
{
	assert(!m_data.empty());
	outseq.clear();
	outseq.reserve(m_data.front().getSequenceLength()
			+ m_data.size() - 1);

	if (m_dir == SENSE) {
		BranchData::const_iterator iter = m_data.begin();
		outseq = iter->decode();
		++iter;
		for (; iter != m_data.end(); ++iter)
			outseq.append(1, iter->getLastBase());
	} else {
		BranchData::const_reverse_iterator iter = m_data.rbegin();
		outseq = iter->decode();
		++iter;
		for (; iter != m_data.rend(); ++iter)
			outseq.append(1, iter->getLastBase());
	}
}

//
// Print the branch
//
void BranchRecord::printBranch(std::ostream& ostr)
{
	for(BranchData::const_iterator iter = m_data.begin(); iter != m_data.end(); ++iter)
	{
		int m = getMultiplicity(*iter);
		ostr << iter->decode() << "," << m << " ";
	}
}

//
// Get the starting iterator
//
BranchDataIter BranchRecord::getStartIter()
{
	return m_data.begin();
}

//
// Get the ending iterator
//
BranchDataIter BranchRecord::getEndIter()
{
	return m_data.end();
}

const PackedSeq& BranchRecord::getSeqByIndex(size_t index) const
{
	assert(index < m_data.size());
	return m_data[index];
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
	PackedSeq first = getFirstSeq();
	PackedSeq last = getLastSeq();
	if (getDirection() == SENSE)
		last.reverseComplement();
	else
		first.reverseComplement();
	return first.compare(last) <= 0;
}
