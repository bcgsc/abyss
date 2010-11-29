#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H 1

#include "ConstString.h"
#include "ContigNode.h"
#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

/** Comparison function used with lower_bound. */
template <typename T>
bool compareLB(const T& a, const typename T::first_type& b)
{
	return a.first < b;
}

/** Comparison function used with upper_bound. */
template <typename T>
bool compareUB(const typename T::first_type& a, const T& b)
{
	return a < b.first;
}

/** A suffix array augmented with a mapped value type. */
class SuffixArray {
  public:
	typedef cstring key_type;
	typedef ContigNode mapped_type;
	typedef std::pair<key_type, mapped_type> value_type;
	typedef std::vector<value_type>::const_iterator const_iterator;

	/** Construct an empty suffix array. */
	SuffixArray(unsigned minOverlap)
		: m_minOverlap(minOverlap), m_dirty(false) { }

	/** Insert the specified sequence into this suffix array. */
	template <typename T>
	void insert(const T& seq, const mapped_type& data)
	{
		m_dirty = true;
		typedef typename T::const_pointer It;
		It last = &seq[seq.size() - m_minOverlap + 1];
		for (It it = &seq[1]; it < last; ++it)
			m_data.push_back(value_type(it, data));
	}

	/** Insert the specified sequence into this suffix array. */
	template <typename T>
	void insert(const std::pair<T, mapped_type>& x)
	{
		insert(x.first, x.second);
	}

	/** Construct the suffix array. */
	void construct()
	{
		if (m_dirty)
			sort(m_data.begin(), m_data.end());
		m_dirty = false;
	}

	/** Find all the elements whose suffix matches the prefix of the
	 * specified query sequence.
	 * @return the range of matches as a pair of iterators
	 */
	template <typename T>
	std::pair<const_iterator, const_iterator> equal_range(
			const T& seq) const
	{
		assert(!m_dirty);
		return std::pair<const_iterator, const_iterator>(
			lower_bound(
				m_data.begin(), m_data.end(), key_type(&seq[0]),
				compareLB<value_type>),
			upper_bound(
				m_data.begin(), m_data.end(), key_type(&seq[0]),
				compareUB<value_type>));
	}

  private:
	unsigned m_minOverlap;
	bool m_dirty;
	std::vector<value_type> m_data;
};

#endif
