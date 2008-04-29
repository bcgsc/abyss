#include "BranchRecord.h"

//
// Constructor
//
BranchRecord::BranchRecord() : m_dir(SENSE), m_state(BS_ACTIVE), m_maxLength(0)
{
	
}

//
// Constructor
//
BranchRecord::BranchRecord(extDirection dir, size_t maxLength) : m_dir(dir), m_state(BS_ACTIVE), m_maxLength(maxLength)
{
	
}

//
// Copy Constructor
//
BranchRecord::BranchRecord(const BranchRecord& other)
{
	m_dir = other.m_dir;
	m_maxLength = other.m_maxLength;
	m_state = other.m_state;
}

//
// Add a single sequence to the branch
//
void BranchRecord::addSequence(const PackedSeq& seq)
{
	m_data.push_back(seq);
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


