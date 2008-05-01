#ifndef BRANCHGROUP_H
#define BRANCHGROUP_H

#include "BranchRecord.h"

enum BranchGroupStatus 
{
	BGS_ACTIVE,
	BGS_JOINED,
	BGS_TOOLONG,
	BGS_LOOPFOUND,
	BGS_TOOMANYBRANCHES
};

typedef std::vector<BranchRecord> BranchGroupData;
class BranchGroup
{
	public:
		BranchGroup(uint64_t id, size_t maxNumBranches);
		
		// Add a branch to the group
		void addBranch(BranchRecord& branch);
		
		// Get the last branch added to the group
		BranchRecord& getLastBranch();
		
		// Get a branch by index
		BranchRecord& getBranch(unsigned int index);
		
		// Get the number of groups
		size_t getNumBranches() const { return m_branches.size(); }
		
		// Check the stop conditions for the branch growth
		BranchGroupStatus checkBubbleStopConditions() const;
		
		// Select a branch to keep (for bubble removal) and return it's index
		size_t selectBranchToKeep() const;
		
		// Iterator accessors 
		
	private:
		BranchGroupData m_branches;
		uint64_t m_id;
		size_t m_maxNumBranches;		
};

#endif
