#ifndef COUNTINGBLOOMFILTERWINDOW_H
#define COUNTINGBLOOMFILTERWINDOW_H 1

#include "Bloom.h"
#include "CascadingBloomFilter.h"
#include <vector>

class CascadingBloomFilterWindow : public CascadingBloomFilter
{
  public:

	/** Constructor.
	 *
	 * @param fullBloomSize size in bits of the containing counting bloom filter
	 * @param startBitPos index of first bit in the window
	 * @param endBitPos index of last bit in the window
	 */
	CascadingBloomFilterWindow(size_t fullBloomSize, size_t startBitPos, size_t endBitPos)
		: m_fullBloomSize(fullBloomSize)
	{
		for (unsigned i = 0; i < MAX_COUNT; i++)
			m_data.push_back(new BloomFilterWindow(m_fullBloomSize, startBitPos, endBitPos));
	}

	/** Add the object with the specified index to this multiset. */
	void insert(size_t index)
	{
		for (unsigned i = 0; i < MAX_COUNT; ++i) {
			assert(m_data.at(i) != NULL);
			if (!(*static_cast<BloomFilterWindow*>(m_data.back()))[index]) {
				static_cast<BloomFilterWindow*>(m_data.back())->insert(index);
				break;
			}
		}
	}

	/** Add the object to this counting multiset. */
	void insert(const Bloom::key_type& key)
	{
		insert(Bloom::hash(key) % m_fullBloomSize);
	}

	void write(std::ostream& out) const
	{
		assert(m_data.back() != NULL);
		out << *static_cast<BloomFilterWindow*>(m_data.back());
	}

	/** Operator for writing the bloom filter to a stream */
	friend std::ostream& operator<<(std::ostream& out, const CascadingBloomFilterWindow& o)
	{
		o.write(out);
		return out;
	}

	/** Get the Bloom filter for a given level */
	BloomFilterWindow& getBloomFilter(unsigned level)
	{
		assert(m_data.at(level) != NULL);
		return *static_cast<BloomFilterWindow*>(m_data.at(level));
	}

  private:
	size_t m_fullBloomSize;
};

#endif
