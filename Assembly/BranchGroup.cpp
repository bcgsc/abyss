#include "BranchGroup.h"

using namespace std;

//
// Constructors
//

BranchGroup::BranchGroup()
	: m_id(0), m_dir(SENSE), m_maxNumBranches(0),
	m_noExt(false), m_status(BGS_ACTIVE), m_branchToKeep(-1)
{
}

BranchGroup::BranchGroup(uint64_t id, extDirection dir,
		size_t maxNumBranches, const PackedSeq &origin)
	: m_id(id), m_dir(dir), m_origin(origin),
	m_maxNumBranches(maxNumBranches), m_noExt(false),
	m_status(BGS_ACTIVE), m_branchToKeep(-1)
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
	m_origin = other.m_origin;
	m_maxNumBranches = other.m_maxNumBranches;
	m_noExt = other.m_noExt;
	m_status = other.m_status;
	m_branchToKeep = other.m_branchToKeep;
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
	m_origin = other.m_origin;
	m_maxNumBranches = other.m_maxNumBranches;
	m_noExt = other.m_noExt;
	m_status = other.m_status;	
	m_branchToKeep = other.m_branchToKeep;
	
	return *this;
}

// Add a branch to the group
BranchRecord& BranchGroup::addBranch(uint64_t id, BranchRecord& branch)
{
	BranchGroupData::iterator newBranch = m_branches.insert(
			pair<uint64_t, BranchRecord>(id, branch)).first;
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
	
	BranchGroupData::const_iterator it = m_branches.begin();
	const PackedSeq& lastSeq = it->second.getLastSeq();
	while (++it != m_branches.end())
		if (it->second.getLastSeq() != lastSeq)
			return m_status = BGS_ACTIVE;
	// All the branches of the bubble have joined.
	selectBranchToKeep();
	return m_status = BGS_JOINED;
}

// 
// Select a branch that survives the bubble removal
//
void BranchGroup::selectBranchToKeep()
{
	assert(m_branchToKeep == -1);
	// Choose the branch with the highest total multiplicity
	// arbitrarily return 0 for now
	int bestMult = 0;
	int bestIndex = -1;
	size_t numBranches = getNumBranches();
	
	for(size_t index = 0; index < numBranches; ++index)
	{
		BranchRecord& branch = getBranch(index);
		int currMult = branch.calculateBranchMultiplicity(true);
		branch.clearMultiplicity();
		if(currMult > bestMult)
		{
			bestMult = currMult;
			bestIndex = index;
		}

		// Remove the last base, which is identical for every branch.
		branch.truncate(branch.getEndIter() - 1);
	}

	assert(bestIndex >= 0);
	m_branchToKeep = bestIndex;
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

/** Return whether this branch is ambiguous at its origin. Also
 * returns false if the origin of the branch has since been deleted.
 */
bool BranchGroup::isAmbiguous(const ISequenceCollection* c) const
{
	if (c->checkFlag(m_origin, SF_DELETE))
		return false;
	ExtensionRecord ext;
	int multiplicity;
	bool found = c->getSeqData(m_origin, ext, multiplicity);
	assert(found);
	(void)found;
	return ext.dir[m_dir].isAmbiguous();
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
