#include "BranchRecord.h"

//
// Constructor
//
BranchRecord::BranchRecord() : m_dir(SENSE), m_state(BS_ACTIVE), m_maxLength(-1), m_loopDetected(false)
{
	
}

//
// Constructor
//
BranchRecord::BranchRecord(extDirection dir, int maxLength) : m_dir(dir), m_state(BS_ACTIVE), m_maxLength(maxLength), m_loopDetected(false)
{
	
}

//
// Copy Constructor
//
BranchRecord::BranchRecord(const BranchRecord& other)
{
	m_data = other.m_data;
	m_set = other.m_set;
	m_dir = other.m_dir;
	m_maxLength = other.m_maxLength;
	m_state = other.m_state;
	m_loopDetected = other.m_loopDetected;
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
	m_set = other.m_set;
	m_dir = other.m_dir;
	m_maxLength = other.m_maxLength;
	m_state = other.m_state;
	m_loopDetected = other.m_loopDetected;	
	return *this;
}

//
// Add a single sequence to the branch
//
void BranchRecord::addSequence(const PackedSeq& seq)
{
	m_data.push_back(seq);
	
	// Detect a loop by checking that the sequence is not already in the branch
	bool unique = m_set.insert(seq).second;
	if(!unique)
	{
		m_loopDetected = true;
	}
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
	return m_set.find(seq) != m_set.end();	
}

//
// If the maxLength == -1, no length check should be performed
//
bool BranchRecord::doLengthCheck() const
{
	return (m_maxLength > -1);
}


//
// Build a contig from a branch
//
void BranchRecord::buildContig(Sequence& outseq) const
{
	// Is there any sequences in the record?
	if(m_data.empty())
	{
		return;
	}
	
	// Calculate the length of the output sequence
	size_t outlen = m_data.front().getSequenceLength() + (m_data.size() - 1);
	// reserve enough room for the entire sequence
	outseq.reserve(outlen);
	
	if(m_dir == SENSE)
	{
		bool first = true;
		for(BranchData::const_iterator iter = m_data.begin(); iter != m_data.end(); ++iter)
		{
			if(first)
			{
				// add the first sequence in its entirety
				outseq.append(iter->decode());
				first = false;
			}
			else
			{
				outseq.append(1, iter->getLastBase());
			}
		}
	}
	else
	{
		// for antisense, start at the end and traverse backwards
		bool first = true;
		for(BranchData::const_reverse_iterator riter = m_data.rbegin(); riter != m_data.rend(); ++riter)
		{
			if(first)
			{
				// add the first sequence in its entirety
				outseq.append(riter->decode());
				first = false;
			}
			else
			{
				outseq.append(1, riter->getLastBase());
			}
		}
	}
	return;
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


