#ifndef COUNTINGBLOOMFILTERWINDOW_H
#define COUNTINGBLOOMFILTERWINDOW_H 1

#include "Bloom.h"
#include "BloomFilter.h"
#include "CascadingBloomFilter.h"
#include <vector>

class CascadingBloomFilterWindow : private CascadingBloomFilter
{
  public:

	/** Constructor.
	 *
	 * @param fullBloomSize size in bits of the containing counting bloom filter
	 * @param startBitPos index of first bit in the window
	 * @param endBitPos index of last bit in the window
	 * @param max_count the maximum count value of the Bloom filter
	 */
	CascadingBloomFilterWindow(size_t fullBloomSize, size_t startBitPos,
			size_t endBitPos, unsigned max_count, size_t hashSeed=0)
		: m_fullBloomSize(fullBloomSize), m_hashSeed(hashSeed)
	{
		m_data.reserve(max_count);
		for (unsigned i = 0; i < max_count; i++)
			m_data.push_back(new BloomFilterWindow(fullBloomSize,
						startBitPos, endBitPos, hashSeed));
	}

	/** Return the size of the bit array. */
	size_t size() const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->size();
	}

	/** Return the number of elements with count >= max_count. */
	size_t popcount() const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->popcount();
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return (double)popcount() / size();
	}

	/** Add the object with the specified index to this multiset. */
	void insert(size_t index)
	{
		for (unsigned i = 0; i < m_data.size(); ++i) {
			assert(m_data.at(i) != NULL);
			if (!(*m_data[i])[index]) {
				m_data[i]->insert(index);
				break;
			}
		}
	}

	/** Add the object to this counting multiset. */
	void insert(const Bloom::key_type& key)
	{
		assert(m_data.back() != NULL);
		insert(Bloom::hash(key, m_hashSeed) % m_fullBloomSize);
	}

	void write(std::ostream& out) const
	{
		assert(m_data.back() != NULL);
		out << *m_data.back();
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
		return *m_data.at(level);
	}

  private:
	size_t m_fullBloomSize;
	size_t m_hashSeed;
	std::vector<BloomFilterWindow*> m_data;
};

#endif
