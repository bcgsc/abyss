#include "BranchGroup.h"

BranchGroup::BranchGroup(uint64_t id, size_t maxNumBranches) : m_id(id), m_maxNumBranches(maxNumBranches)
{
	
	
}

// Add a branch to the group
void BranchGroup::addBranch(BranchRecord& branch)
{
	m_branches.push_back(branch);
}

// Get the last branch added to the group
BranchRecord& BranchGroup::getLastBranch()
{
	assert(!m_branches.empty());
	return m_branches.back();
}

// Get a branch by index
BranchRecord& BranchGroup::getBranch(unsigned int index)
{
	assert(index < m_branches.size());
	return m_branches[index];
}

//		
// Check the stop conditions for the bubble growth
//
BranchGroupStatus BranchGroup::checkBubbleStopConditions() const
{
	// Check if there are too many branches
	if(m_branches.size() > m_maxNumBranches)
	{
		return BGS_TOOMANYBRANCHES;
	}
	
	// Check if any branches are too long or any sequence has a loop
	for(BranchGroupData::const_iterator iter = m_branches.begin(); iter != m_branches.end(); ++iter)
	{
		if(iter->getLength() > iter->getMaxLength())
		{
			return BGS_TOOLONG;
		}
		
		if(iter->hasLoop())
		{
			return BGS_LOOPFOUND;
		}
	}
	
	// Check if the sequences have joined back together
	
	// For each sequence check if the last added sequence is present in ALL the other branches
	// This indidates the branches have joined back together
	// It is not necessary that the bubbles are the same length
	for(BranchGroupData::const_iterator iter = m_branches.begin(); iter != m_branches.end(); ++iter)
	{
		const PackedSeq& lastSeq = iter->getLastSeq();
		bool allFound = true;
		for(BranchGroupData::const_iterator otherIter = m_branches.begin(); otherIter != m_branches.end(); ++otherIter)
		{
			// Skip selfcheck
			if(iter == otherIter)
			{
				continue;
			}
			
			// Check if the sequence is in this branch
			if(!otherIter->exists(lastSeq))
			{
				allFound = false;
				break;
			}
		}
		
		if(allFound)
		{
			return BGS_JOINED;
		}
	}
	
	// If we reached here, no stop condition was met and no join was found so the extension can continue
	return BGS_ACTIVE;
}

// 
// Select a branch that survives the bubble removal
//
size_t BranchGroup::selectBranchToKeep() const
{
	// arbitrarily return 0 for now
	return 0;
}
