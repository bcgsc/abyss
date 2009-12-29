#include "BranchGroup.h"

using namespace std;

// Check the stop conditions for the bubble growth
BranchGroupStatus BranchGroup::updateStatus()
{
	assert(m_branches.size() <= m_maxNumBranches);

	if (m_status != BGS_ACTIVE)
		return m_status;

	// Check if the no extension flag is set
	if(m_noExt)
	{
		m_status = BGS_NOEXT;
		return m_status;
	}

	// Check if any branches are too long or any sequence has a loop
	for(BranchGroupData::const_iterator iter = m_branches.begin(); iter != m_branches.end(); ++iter)
	{
		if (iter->getLength() > iter->getMaxLength()) {
			m_status = BGS_TOOLONG;
			return m_status;
		}

		if (iter->hasLoop()) {
			m_status = BGS_LOOPFOUND;
			return m_status;
		}
	}

	BranchGroupData::const_iterator it = m_branches.begin();
	const Kmer& lastSeq = it->getLastSeq();
	while (++it != m_branches.end())
		if (it->getLastSeq() != lastSeq)
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
		branch.truncate(branch.end() - 1);
	}

	assert(bestIndex >= 0);
	m_branchToKeep = bestIndex;
}

/** Return whether any branches of this group are active. */
bool BranchGroup::isActive() const
{
	for (BranchGroupData::const_iterator it = m_branches.begin();
			it != m_branches.end(); ++it)
		if (it->isActive())
			return true;
	return false;
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
	// Get fresh data from the collection to check that this bubble
	// does in fact still exist.
	const KmerData& data = c->getSeqAndData(m_origin).second;
	return data.deleted() ? false : data.isAmbiguous(m_dir);
}
