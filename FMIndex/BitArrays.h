#ifndef BITARRAYS_H
#define BITARRAYS_H 1

#include "bit_array.h"
#include <algorithm>
#include <cassert>
#include <cstdlib> // for abort
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
  public:

/** Count the occurences of the symbols of s. */
template <typename T>
void Init(const std::vector<T>& s)
{
	assert(!s.empty());
	m_data.clear();
	unsigned n = *std::max_element(s.begin(), s.end()) + 1;
	assert(n > 0);
	assert(n < std::numeric_limits<T>::max());
	m_data.resize(n, wat_array::BitArray(s.size()));

	typedef typename std::vector<T>::const_iterator It;
	size_t i = 0;
	for (It it = s.begin(); it != s.end(); ++it, ++i) {
		T c = *it;
		assert(c < m_data.size());
		m_data[c].SetBit(1, i);
	}
}

/** Return the size of the string. */
size_t size() const
{
	assert(!m_data.empty());
	return m_data.front().length();
}

/** Return the size of the string. */
size_t length() const
{
	return size();
}

/** Return the count of symbol c in s[0, i). */
size_t Rank(uint8_t c, size_t i) const
{
	return m_data[c].Rank(1, i);
}

/** Return s[i]. */
uint8_t Lookup(size_t i) const
{
	for (Data::const_iterator it = m_data.begin();
			it != m_data.end(); ++it)
		if (it->Lookup(i))
			return it - m_data.begin();
	assert(false);
	abort();
}

/** Load this data structure. */
void Load(std::istream& in)
{
	typedef uint8_t T;
	m_data.clear();
	uint32_t n = 0;
	if (!in.read(reinterpret_cast<char*>(&n), sizeof n))
		return;
	assert(n > 0);
	assert(n < std::numeric_limits<T>::max());
	m_data.resize(n);
	for (Data::iterator it = m_data.begin();
			it != m_data.end(); ++it)
		it->Load(in);
}

/** Store this data structure. */
void Save(std::ostream& out) const
{
	uint32_t n = m_data.size();
	out.write(reinterpret_cast<char*>(&n), sizeof n);
	for (Data::const_iterator it = m_data.begin();
			it != m_data.end(); ++it)
		it->Save(out);
}

  private:
	typedef std::vector<wat_array::BitArray> Data;
	Data m_data;
};

#endif
