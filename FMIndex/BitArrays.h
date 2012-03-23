#ifndef BITARRAYS_H
#define BITARRAYS_H 1

#include "bit_array.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <istream>
#include <limits> // for numeric_limits
#include <ostream>
#include <stdint.h>
#include <vector>

/** Store a string of symbols from a small alphabet using a vector of
 * BitArray.
 */
class BitArrays
{
	/** A symbol. */
	typedef uint8_t T;

	/** The sentinel symbol. */
	static T SENTINEL() { return std::numeric_limits<T>::max(); }

  public:

/** Count the occurrences of the symbols of [first, last). */
template<typename It>
void assign(It first, It last)
{
	assert(first < last);
	m_data.clear();

	// Determine the size of the alphabet ignoring the sentinel.
	T n = 0;
	for (It it = first; it != last; ++it)
		if (*it != SENTINEL())
			n = std::max(n, *it);
	n++;

	assert(n < std::numeric_limits<T>::max());
	m_data.resize(n, wat_array::BitArray(last - first));

	size_t i = 0;
	for (It it = first; it != last; ++it, ++i) {
		T c = *it;
		if (c == SENTINEL())
			continue;
		assert(c < m_data.size());
		m_data[c].SetBit(1, i);
	}

	std::for_each(m_data.begin(), m_data.end(),
			std::mem_fun_ref(&wat_array::BitArray::Build));
}

/** Return the size of the string. */
size_t size() const
{
	assert(!m_data.empty());
	return m_data.front().length();
}

/** Return the number of occurrences of the specified symbol. */
size_t count(T c) const
{
	return m_data[c].one_num();
}

/** Return the count of symbol c in s[0, i). */
size_t rank(T c, size_t i) const
{
	return m_data[c].Rank(1, i);
}

/** Return the symbol at the specified position. */
T at(size_t i) const
{
	assert(!m_data.empty());
	assert(i < m_data.front().length());
	for (Data::const_iterator it = m_data.begin();
			it != m_data.end(); ++it)
		if (it->Lookup(i))
			return it - m_data.begin();
	return std::numeric_limits<T>::max();
}

/** Store this data structure. */
friend std::ostream& operator<<(std::ostream& out, const BitArrays& o)
{
	uint32_t n = o.m_data.size();
	out.write(reinterpret_cast<char*>(&n), sizeof n);
	for (Data::const_iterator it = o.m_data.begin();
			it != o.m_data.end(); ++it)
		it->Save(out);
	return out;
}

/** Load this data structure. */
friend std::istream& operator>>(std::istream& in, BitArrays& o)
{
	o.m_data.clear();
	uint32_t n = 0;
	if (!in.read(reinterpret_cast<char*>(&n), sizeof n))
		return in;
	assert(n > 0);
	assert(n < std::numeric_limits<T>::max());
	o.m_data.resize(n);
	for (Data::iterator it = o.m_data.begin();
			it != o.m_data.end(); ++it)
		it->Load(in);
	return in;
}

  private:
	typedef std::vector<wat_array::BitArray> Data;
	Data m_data;
};

#endif
