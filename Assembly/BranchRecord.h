#ifndef BRANCHRECORD_H
#define BRANCHRECORD_H 1

#include "Kmer.h"
#include "KmerData.h"
#include <utility>
#include <vector>

enum BranchState
{
	BS_ACTIVE, // can be extended
	BS_NOEXT, // the branch has ended because of a lack of sequence to extend to
	BS_AMBI_SAME, // the branch has ended because the extension from this branch is ambigious
	BS_AMBI_OPP, // the branch has ended because the extension to this branch is ambigiuous
	BS_TOO_LONG, // the branch is too long
};

class BranchRecord
{
	public:
		typedef std::pair<Kmer, KmerData> value_type;
		typedef std::vector<value_type> BranchData;
		typedef BranchData::iterator iterator;
		typedef BranchData::const_iterator const_iterator;

		BranchRecord()
			: m_dir(SENSE), m_state(BS_ACTIVE), m_maxLength(-1) { }
		BranchRecord(extDirection dir, int maxLength)
			: m_dir(dir), m_state(BS_ACTIVE), m_maxLength(maxLength)
			{ }

		operator Sequence() const;

		/** Return true if this sequence has no elements. */
		bool empty() const { return m_data.empty(); }

		/** Return the number of elements. */
		size_t size() const { return m_data.size(); }

		/** Add the element x at the end. */
		void push_back(const value_type& x) { m_data.push_back(x); }

		/** Remove the last k-mer. */
		void pop_back()
		{
			assert(!m_data.empty());
			m_data.pop_back();
		}

		/** Return the first element. */
		const value_type& front() const
		{
			assert(!m_data.empty());
			return m_data.front();
		}

		/** Return the last element. */
		const value_type& back() const
		{
			assert(!m_data.empty());
			return m_data.back();
		}

		/** Terminate this branch with the specified reason. */
		void terminate(BranchState reason)
		{
			assert(reason != BS_ACTIVE);
			m_state = reason;
		}

		/** Return whether this branch is active. */
		bool isActive() const { return m_state == BS_ACTIVE; }

		/** Return the state of this branch. */
		BranchState getState() const { return m_state;	}

		/** Return the direction of this branch. */
		extDirection getDirection() const { return m_dir; }

		/** Set the properties of the last element. */
		void setData(const value_type& o)
		{
			assert(m_data.back().first == o.first);
			m_data.back().second = o.second;
		}

		iterator begin() { return m_data.begin(); }
		iterator end() { return m_data.end(); }
		const_iterator begin() const { return m_data.begin(); }
		const_iterator end() const { return m_data.end(); }

		// Return the maximum branch length
		size_t getMaxLength() const { return m_maxLength; }

		// check if a sequence exists in the branch record
		bool exists(const Kmer& seq) const;

		/** Return whether the branch is too long. */
		bool isTooLong() const;

		int calculateBranchMultiplicity() const;

		bool isCanonical() const;

	private:
		BranchData m_data;
		extDirection m_dir;
		BranchState m_state;
		int m_maxLength;
};

#endif
