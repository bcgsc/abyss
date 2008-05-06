#include "BranchGroup.h"


//
// Constructors
//

BranchGroup::BranchGroup() : m_id(0), m_dir(SENSE), m_maxNumBranches(0), m_noExt(false), m_status(BGS_ACTIVE)
{
	
	
}

BranchGroup::BranchGroup(uint64_t id, extDirection dir, size_t maxNumBranches) : m_id(id), m_dir(dir), m_maxNumBranches(maxNumBranches), m_noExt(false), m_status(BGS_ACTIVE)
{
	
	
}

//
// Copy constructor
//
BranchGroup::BranchGroup(const BranchGroup& other)
{
	m_branches = other.m_branches;
	m_id = other.m_id;
	m_dir = other.m_dir;
	m_maxNumBranches = other.m_maxNumBranches;
	m_noExt = other.m_noExt;
	m_status = other.m_status;
	
}

//
// Assignment operator
//
BranchGroup& BranchGroup::operator=(const BranchGroup& other)
{
	// Detect self assignment
	if (this == &other)
	{
		return *this;
	}
			
	m_branches = other.m_branches;
	m_id = other.m_id;
	m_dir = other.m_dir;
	m_maxNumBranches = other.m_maxNumBranches;
	m_noExt = other.m_noExt;
	m_status = other.m_status;	
	
	return *this;
	
}

// Add a branch to the group
BranchRecord& BranchGroup::addBranch(uint64_t id, BranchRecord& branch)
{
	BranchGroupData::iterator newBranch = m_branches.insert(std::pair<uint64_t, BranchRecord>(id, branch)).first;
	return newBranch->second;
	
}

//
// Get a branch by ID
//
BranchRecord& BranchGroup::getBranch(uint64_t id)
{
	BranchGroupData::iterator branch = m_branches.find(id);
	// should never fail
	assert(branch != m_branches.end());
	return branch->second;
}

//		
// Check the stop conditions for the bubble growth
//
BranchGroupStatus BranchGroup::updateStatus()
{
	// Check if the no extension flag is set
	if(m_noExt)
	{
		m_status = BGS_NOEXT;
		return m_status;
	}
	
	// Check if there are too many branches
	if(m_branches.size() > m_maxNumBranches)
	{
		m_status = BGS_TOOMANYBRANCHES;
		return m_status;
	}
	
	// Check if any branches are too long or any sequence has a loop
	for(BranchGroupData::const_iterator iter = m_branches.begin(); iter != m_branches.end(); ++iter)
	{
		if(iter->second.getLength() > iter->second.getMaxLength())
		{
			m_status = BGS_TOOLONG;
			return m_status;
		}
		
		if(iter->second.hasLoop())
		{
			m_status = BGS_LOOPFOUND;
			return m_status;
		}
	}
	
	// Check if the sequences have joined back together
	
	// For each sequence check if the last added sequence is present in ALL the other branches
	// This indidates the branches have joined back together
	// It is not necessary that the bubbles are the same length
	for(BranchGroupData::const_iterator iter = m_branches.begin(); iter != m_branches.end(); ++iter)
	{
		const PackedSeq& lastSeq = iter->second.getLastSeq();
		bool allFound = true;
		for(BranchGroupData::const_iterator otherIter = m_branches.begin(); otherIter != m_branches.end(); ++otherIter)
		{
			// Skip selfcheck
			if(iter == otherIter)
			{
				continue;
			}
			
			// Check if the sequence is in this branch
			if(!otherIter->second.exists(lastSeq))
			{
				allFound = false;
				break;
			}
		}
		
		if(allFound)
		{
			m_status = BGS_JOINED;
			return m_status;
		}
	}
	
	// If we reached here, no stop condition was met and no join was found so the extension can continue
	m_status = BGS_ACTIVE;
	return m_status;
}

// 
// Select a branch that survives the bubble removal
//
size_t BranchGroup::selectBranchToKeep() const
{
	// arbitrarily return 0 for now
	return 0;
}

//
// Is the branch group active
//
bool BranchGroup::isActive() const
{
	// Are any of the branches active?
	bool active = false;	
	for(BranchGroupData::const_iterator iter = m_branches.begin(); iter != m_branches.end(); ++iter)
	{
		active = (active || iter->second.isActive());	
	}
	return active;
}

//
// Is the branch extendable?
//
bool BranchGroup::isExtendable()
{
	// A group is extendable when all the branches are the same length (they are lockstepped for growth) and the noextension flag is not set
	bool sameLength = true;
	size_t numBranches = getNumBranches();
	for(size_t index = 0; index < numBranches - 1; index++)
	{
		sameLength = (sameLength && (getBranch(index).getLength() == getBranch(index+1).getLength()));
	}
	
	// If the branches are all the same length and the noext flag is not set, return true
	bool ret = sameLength && !isNoExt();
	
	return ret;
}

//
// Print the sizes of the branches
//
void BranchGroup::printBranches() const
{
	for(BranchGroupData::const_iterator iter = m_branches.begin(); iter != m_branches.end(); ++iter)
	{
		printf("	%zu\n", iter->second.getLength());
	}	
	
}
