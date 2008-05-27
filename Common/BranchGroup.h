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

typedef std::map<uint64_t, BranchRecord> BranchGroupData;

class BranchGroup
{
	public:
		BranchGroup();
		BranchGroup(uint64_t id, extDirection dir, size_t maxNumBranches);
		BranchGroup(const BranchGroup& other);
		
		BranchGroup& operator=(const BranchGroup& other);
		
		// Add a branch to the group
		BranchRecord& addBranch(uint64_t id, BranchRecord& branch);
		
		// Get the branch corresponding to an id
		BranchRecord& getBranch(uint64_t id);
		
		// Get the number of groups
		size_t getNumBranches() const { return m_branches.size(); }
		
		// Set the maximum number of branches.
		void setMaxNumBranches(size_t maxNumBranches)
		{
			m_maxNumBranches = maxNumBranches;
		}

		// Check the stop conditions for the branch growth
		BranchGroupStatus updateStatus();
		
		// return the current status of the branch
		BranchGroupStatus getStatus() const { return m_status; }
		
		// Select a branch to keep (for bubble removal) and return it's index
		size_t selectBranchToKeep();
		
		// set the no extension flag
		void setNoExtension() { m_noExt = true; }
		
		// check if the group is active
		bool isActive() const;
		
		// is the no extension flag set?
		bool isNoExt() const { return m_noExt; }
		
		// check if the group is ready for another round of extension (all the branches are the same length)
		bool isExtendable();
		
		// return the direction of growth
		extDirection getDirection() const { return m_dir; }
		
		// print the branches that make up this group
		void printBranches() const;
		
		// Iterator accessors 
		BranchGroupData::iterator getStartIter() { return m_branches.begin(); }
		BranchGroupData::iterator getEndIter() { return m_branches.end(); }
		
	private:
		BranchGroupData m_branches;
		uint64_t m_id;
		extDirection m_dir;
		size_t m_maxNumBranches;
		bool m_noExt;
		BranchGroupStatus m_status;
};

#endif
