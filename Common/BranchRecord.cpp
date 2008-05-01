#include "BranchRecord.h"

//
// Constructor
//
BranchRecord::BranchRecord() : m_dir(SENSE), m_state(BS_ACTIVE), m_maxLength(0), m_loopDetected(false)
{
	
}

//
// Constructor
//
BranchRecord::BranchRecord(extDirection dir, size_t maxLength) : m_dir(dir), m_state(BS_ACTIVE), m_maxLength(maxLength), m_loopDetected(false)
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


