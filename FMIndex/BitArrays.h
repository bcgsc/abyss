#ifndef BITARRAYS_H
#define BITARRAYS_H 1

#include "bit_array.h"
#include <algorithm>
#include <cassert>
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

/** Count the occurences of the symbols of s. */
void assign(const std::vector<T>& s)
{
	assert(!s.empty());
	m_data.clear();

	// Determine the size of the alphabet ignoring the sentinel.
	T n = 0;
	for (std::vector<T>::const_iterator it = s.begin();
			it != s.end(); ++it)
		if (*it != SENTINEL())
			n = std::max(n, *it);
	n++;

	assert(n < std::numeric_limits<T>::max());
	m_data.resize(n, wat_array::BitArray(s.size()));

	typedef std::vector<T>::const_iterator It;
	size_t i = 0;
	for (It it = s.begin(); it != s.end(); ++it, ++i) {
		T c = *it;
		if (c == SENTINEL())
			continue;
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

/** Load this data structure. */
void load(std::istream& in)
{
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
void save(std::ostream& out) const
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
