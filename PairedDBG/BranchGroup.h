#ifndef PAIREDDBG_BRANCHGROUP_H
#define PAIREDDBG_BRANCHGROUP_H 1

#include "BranchRecord.h"
#include "SequenceCollection.h"
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
			: m_dir(SENSE), m_maxNumBranches(0),
			m_noExt(false), m_status(BGS_ACTIVE)
			{ }

		BranchGroup(extDirection dir, size_t maxNumBranches,
				const Kmer &origin)
			: m_dir(dir), m_origin(origin),
			m_maxNumBranches(maxNumBranches), m_noExt(false),
			m_status(BGS_ACTIVE)
		{
			m_branches.reserve(m_maxNumBranches);
		}

		BranchGroup(extDirection dir, size_t maxNumBranches,
				const Kmer &origin, const BranchRecord& branch)
			: m_dir(dir), m_origin(origin),
			m_maxNumBranches(maxNumBranches), m_noExt(false),
			m_status(BGS_ACTIVE)
		{
			m_branches.reserve(m_maxNumBranches);
			m_branches.push_back(branch);
		}

		BranchGroup(const BranchGroup& o)
			: m_branches(o.m_branches), m_dir(o.m_dir),
			m_origin(o.m_origin),
			m_maxNumBranches(o.m_maxNumBranches),
			m_noExt(o.m_noExt), m_status(o.m_status)
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
		void addBranch(const BranchRecord& branch,
				const Kmer& kmer)
		{
			if (m_branches.size() < m_maxNumBranches)
				addBranch(branch).push_back(
						std::make_pair(kmer, KmerData()));
			else
				m_status = BGS_TOOMANYBRANCHES;
		}

		/** Return the specified branch. */
		BranchRecord& operator [](unsigned id)
		{
			return m_branches[id];
		}

		/** Return the number of branches in this group. */
		size_t size() const { return m_branches.size(); }

		/** Return whether a branch contains the specified k-mer at
		 * the index i. */
		bool exists(unsigned i, const Kmer& kmer) const
		{
			for (BranchGroupData::const_iterator it
					= m_branches.begin();
					it != m_branches.end(); ++it)
				if (it->exists(i, kmer))
					return true;
			return false;
		}

		BranchGroupStatus updateStatus(unsigned maxLength);

		// return the current status of the branch
		BranchGroupStatus getStatus() const { return m_status; }

		// set the no extension flag
		void setNoExtension() { m_noExt = true; }

		bool isActive() const;

		// is the no extension flag set?
		bool isNoExt() const { return m_noExt; }

		bool isExtendable();

		// return the direction of growth
		extDirection getDirection() const { return m_dir; }

		iterator begin() { return m_branches.begin(); }
		iterator end() { return m_branches.end(); }
		const_iterator begin() const { return m_branches.begin(); }
		const_iterator end() const { return m_branches.end(); }

		bool isAmbiguous(const SequenceCollectionHash& c) const;

	private:
		BranchGroup& operator =(const BranchGroup& o);

		BranchGroupData m_branches;
		extDirection m_dir;
 		Kmer m_origin;
		size_t m_maxNumBranches;
		bool m_noExt;
		BranchGroupStatus m_status;
};

namespace std {
	template <>
	inline void swap(BranchGroup&, BranchGroup&) { assert(false); }
}

#endif
