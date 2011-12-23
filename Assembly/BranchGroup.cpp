#include "BranchGroup.h"
#include "Algorithms.h"
#include <algorithm>
#include <functional>

using namespace std;

// Check the stop conditions for the bubble growth
BranchGroupStatus BranchGroup::updateStatus(unsigned maxLength)
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
	for (BranchGroupData::const_iterator iter = m_branches.begin();
			iter != m_branches.end(); ++iter) {
		if (iter->isTooLong(maxLength)) {
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
	// Remove the last base, which is identical for every branch.
	for_each(m_branches.begin(), m_branches.end(),
			mem_fun_ref(&BranchRecord::pop_back));

	// Sort the branches by coverage.
	sort_by_transform(m_branches.begin(), m_branches.end(),
			mem_fun_ref(&BranchRecord::calculateBranchMultiplicity));
	reverse(m_branches.begin(), m_branches.end());

	return m_status = BGS_JOINED;
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
bool BranchGroup::isAmbiguous(const SequenceCollectionHash& g) const
{
	// Get fresh data from the collection to check that this bubble
	// does in fact still exist.
	const KmerData& data = g.getSeqAndData(m_origin).second;
	return data.deleted() ? false : data.isAmbiguous(m_dir);
}
