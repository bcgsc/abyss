#ifndef BRANCHGROUP_H
#define BRANCHGROUP_H

#include "BranchRecord.h"

enum BranchGroupStatus 
{
	BGS_ACTIVE,
	BGS_NOEXT,
	BGS_JOINED,
	BGS_TOOLONG,
	BGS_LOOPFOUND,
	BGS_TOOMANYBRANCHES
};

typedef std::vector<BranchRecord> BranchGroupData;
class BranchGroup
{
	public:
		BranchGroup(uint64_t id, extDirection dir, size_t maxNumBranches);
		
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
		
		// set the no extension flag
		void setNoExtension() { m_noExt = true; }
		
		extDirection getDirection() { return m_dir; }
		// Iterator accessors 
		
	private:
		BranchGroupData m_branches;
		uint64_t m_id;
		extDirection m_dir;
		size_t m_maxNumBranches;
		bool m_noExt;
};

#endif
