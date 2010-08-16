#include "BranchGroup.h"
#include <algorithm>

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
		if (iter->size() > iter->getMaxLength()) {
			m_status = BGS_TOOLONG;
			return m_status;
		}
	}

	BranchGroupData::const_iterator it = m_branches.begin();
	const Kmer& lastSeq = it->back().first;
	while (++it != m_branches.end())
		if (it->back().first != lastSeq)
			return m_status = BGS_ACTIVE;

	// All the branches of the bubble have joined.
	sortByCoverage();
	return m_status = BGS_JOINED;
}

/** Return whether branch a has higher coverage than branch b. */
static bool compareCoverage(
		const BranchRecord& a, const BranchRecord& b)
{
	return a.getBranchMultiplicity() > b.getBranchMultiplicity();
}

/** Sort the branches by coverage. */
void BranchGroup::sortByCoverage()
{
	// Sum up the coverage for each branch.
	for (BranchGroupData::iterator it = m_branches.begin();
			it != m_branches.end(); ++it) {
		// Remove the last base, which is identical for every branch.
		it->pop_back();
		it->calculateBranchMultiplicity();
	}

	// Sort by coverage.
	sort(m_branches.begin(), m_branches.end(), compareCoverage);
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

/** Return whether this branch is extendable. */
bool BranchGroup::isExtendable()
{
	if (m_noExt)
		return false;

	// A group is extendable when all the branches are the same
	// length. All the branches are lockstepped for growth.
	BranchGroupData::iterator it = m_branches.begin();
	unsigned length = it++->size();
	for (; it != m_branches.end(); ++it)
		if (it->size() != length)
			return false;
	return true;
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
