#ifndef BRANCHGROUP_H
#define BRANCHGROUP_H 1

#include "Common/Algorithms.h"
#include "Common/Exception.h"
#include <algorithm> // for swap
#include <map>
#include <utility>

enum BranchGroupStatus
{
	BGS_ACTIVE,
	BGS_NOEXT,
	BGS_JOINED,
	BGS_TOOLONG,
	BGS_TOOMANYBRANCHES
};

/** A container of BranchRecord. */
class BranchGroup
{
  public:
	typedef std::vector<BranchRecord> BranchGroupData;
	typedef BranchGroupData::iterator iterator;
	typedef BranchGroupData::const_iterator const_iterator;

	BranchGroup()
	  : m_dir(SENSE)
	  , m_maxNumBranches(0)
	  , m_noExt(false)
	  , m_status(BGS_ACTIVE)
	{}

	BranchGroup(extDirection dir, size_t maxNumBranches, const BranchRecord::V& origin)
	  : m_dir(dir)
	  , m_origin(origin)
	  , m_maxNumBranches(maxNumBranches)
	  , m_noExt(false)
	  , m_status(BGS_ACTIVE)
	{
		m_branches.reserve(m_maxNumBranches);
	}

	BranchGroup(
	    extDirection dir,
	    size_t maxNumBranches,
	    const BranchRecord::V& origin,
	    const BranchRecord& branch)
	  : m_dir(dir)
	  , m_origin(origin)
	  , m_maxNumBranches(maxNumBranches)
	  , m_noExt(false)
	  , m_status(BGS_ACTIVE)
	{
		m_branches.reserve(m_maxNumBranches);
		m_branches.push_back(branch);
	}

	BranchGroup(const BranchGroup& o)
	  : m_branches(o.m_branches)
	  , m_dir(o.m_dir)
	  , m_origin(o.m_origin)
	  , m_maxNumBranches(o.m_maxNumBranches)
	  , m_noExt(o.m_noExt)
	  , m_status(o.m_status)
	{
		m_branches.reserve(m_maxNumBranches);
	}

	/** Add a branch to this group. */
	BranchRecord& addBranch(const BranchRecord& branch)
	{
		assert(m_branches.size() < m_maxNumBranches);
		m_branches.push_back(branch);
		return m_branches.back();
	}

	/** Add a branch to this group and extend the new branch with
	 * the given k-mer. */
	void addBranch(const BranchRecord& branch, const BranchRecord::V& kmer)
	{
		if (m_branches.size() < m_maxNumBranches)
			addBranch(branch).push_back(std::make_pair(kmer, BranchRecord::VP()));
		else
			m_status = BGS_TOOMANYBRANCHES;
	}

	/** Return the specified branch. */
	BranchRecord& operator[](unsigned id) { return m_branches[id]; }

	/** Return the number of branches in this group. */
	size_t size() const { return m_branches.size(); }

	/** Return whether a branch contains the specified k-mer at
	 * the index i. */
	bool exists(unsigned i, const BranchRecord::V& kmer) const
	{
		for (BranchGroupData::const_iterator it = m_branches.begin(); it != m_branches.end(); ++it)
			if (it->exists(i, kmer))
				return true;
		return false;
	}

	// return the current status of the branch
	BranchGroupStatus getStatus() const { return m_status; }

	// set the no extension flag
	void setNoExtension() { m_noExt = true; }

	// is the no extension flag set?
	bool isNoExt() const { return m_noExt; }

	// return the direction of growth
	extDirection getDirection() const { return m_dir; }

	iterator begin() { return m_branches.begin(); }
	iterator end() { return m_branches.end(); }
	const_iterator begin() const { return m_branches.begin(); }
	const_iterator end() const { return m_branches.end(); }

	// Check the stop conditions for the bubble growth
	BranchGroupStatus updateStatus(unsigned maxLength)
	{
		assert(m_branches.size() <= m_maxNumBranches);

		if (m_status != BGS_ACTIVE)
			return m_status;

		// Check if the no extension flag is set
		if (m_noExt) {
			m_status = BGS_NOEXT;
			return m_status;
		}

		// Check if any branches are too long or any sequence has a loop
		for (BranchGroupData::const_iterator iter = m_branches.begin(); iter != m_branches.end();
		     ++iter) {
			if (iter->isTooLong(maxLength)) {
				m_status = BGS_TOOLONG;
				return m_status;
			}
		}

		BranchGroupData::const_iterator it = m_branches.begin();
		const BranchRecord::V& lastSeq = it->back().first;
		while (++it != m_branches.end())
			if (it->back().first != lastSeq)
				return m_status = BGS_ACTIVE;

		// All the branches of the bubble have joined.
		// Remove the last base, which is identical for every branch.
		std::for_each(
		    m_branches.begin(), m_branches.end(), [](BranchRecord& b) { return b.pop_back(); });

		// Sort the branches by coverage.
		std::function<int(const BranchRecord&)> lambda = [](const BranchRecord& b) {
			return b.calculateBranchMultiplicity();
		};
		sort_by_transform(m_branches.begin(), m_branches.end(), lambda);
		reverse(m_branches.begin(), m_branches.end());

		return m_status = BGS_JOINED;
	}

	/** Return whether any branches of this group are active. */
	bool isActive() const
	{
		for (BranchGroupData::const_iterator it = m_branches.begin(); it != m_branches.end(); ++it)
			if (it->isActive())
				return true;
		return false;
	}

	/** Return whether this branch is extendable. */
	bool isExtendable()
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
	bool isAmbiguous(const SequenceCollectionHash& g) const
	{
		// Get fresh data from the collection to check that this bubble
		// does in fact still exist.
		const BranchRecord::VP& data = g.getSeqAndData(m_origin).second;
		return data.deleted() ? false : data.isAmbiguous(m_dir);
	}

  private:
	BranchGroup& operator=(const BranchGroup& o);

	BranchGroupData m_branches;
	extDirection m_dir;
	BranchRecord::V m_origin;
	size_t m_maxNumBranches;
	bool m_noExt;
	BranchGroupStatus m_status;
};

#endif
