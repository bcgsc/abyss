#ifndef BRANCHRECORD_H
#define BRANCHRECORD_H 1

#include "Kmer.h"
#include "KmerData.h"
#include <algorithm> // for swap
#include <cassert>
#include <utility>
#include <vector>

enum BranchState
{
	// The branch can be extended.
	BS_ACTIVE,
	// The branch has ended because of a lack of sequence to extend to
	BS_NOEXT,
	// The branch has ended because the extension from this branch is
	// ambigious.
	BS_AMBI_SAME,
	// The branch has ended because the extension to this branch is
	// ambigiuous.
	BS_AMBI_OPP,
	// The branch is too long.
	BS_TOO_LONG,
};

/** A sequence of Kmer. */
class BranchRecord
{
	public:
		typedef Kmer V;
		typedef KmerData VP;
		typedef std::pair<V, VP> value_type;
		typedef std::vector<value_type> BranchData;
		typedef BranchData::iterator iterator;
		typedef BranchData::const_iterator const_iterator;

		BranchRecord() : m_dir(SENSE), m_state(BS_ACTIVE) { }

		explicit BranchRecord(extDirection dir)
			: m_dir(dir), m_state(BS_ACTIVE) { }

		void swap(BranchRecord& o)
		{
			std::swap(m_data, o.m_data);
			std::swap(m_dir, o.m_dir);
			std::swap(m_state, o.m_state);
		}

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

		/** Return true if the k-mer at position i is the specified
		 * k-mer. */
		bool exists(unsigned i, const V& kmer) const
		{
			assert(i < m_data.size());
			return m_data[i].first == kmer;
		}

		/** Return true if this branch is longer than maxLength. */
		bool isTooLong(unsigned maxLength) const
		{
			return size() > maxLength;
		}

/** Calculate the total multiplicity of this branch. */
int calculateBranchMultiplicity() const
{
	assert(!m_data.empty());
	int total = 0;
	for (BranchData::const_iterator it = m_data.begin();
			it != m_data.end(); ++it) {
		int m = it->second.getMultiplicity();
		assert(m >= 0);
		total += m;
	}
	assert(total >= 0);
	return total;
}

/**
 * Return whether this branch is the canonical representation of the
 * contig that it represents. A contig has two ends, and the contig
 * is output starting from the lexicographically smaller end.
 */
bool isCanonical() const
{
	assert(size() > 1);
	Kmer first = front().first;
	Kmer last = back().first;
	if (getDirection() == SENSE)
		last.reverseComplement();
	else
		first.reverseComplement();
	assert(first != last);
	return first < last;
}

/** Build a contig from a branch. */
operator Sequence() const
{
	assert(!m_data.empty());
	Sequence outseq;
	outseq.reserve(m_data.front().first.length() + m_data.size() - 1);

	if (m_dir == SENSE) {
		BranchData::const_iterator iter = m_data.begin();
		outseq = iter->first.str();
		++iter;
		for (; iter != m_data.end(); ++iter)
			outseq.append(1, iter->first.getLastBaseChar());
	} else {
		BranchData::const_reverse_iterator iter = m_data.rbegin();
		outseq = iter->first.str();
		++iter;
		for (; iter != m_data.rend(); ++iter)
			outseq.append(1, iter->first.getLastBaseChar());
	}
	return outseq;
}

	private:
		BranchData m_data;
		extDirection m_dir;
		BranchState m_state;
};

namespace std {
	template <>
	inline void swap(BranchRecord& a, BranchRecord& b)
	{
		a.swap(b);
	}
}

#endif
