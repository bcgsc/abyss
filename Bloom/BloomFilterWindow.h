#ifndef BLOOMFILTERWINDOW_H
#define BLOOMFILTERWINDOW_H 1

#include "Bloom.h"
#include "BloomFilter.h"
#include "Common/HashFunction.h"
#include "Common/Kmer.h"
#include "Common/IOUtil.h"
#include <algorithm>
#include <vector>
#include <iostream>

//TODO many copied functions from bloomfilter.h despite inheritance need to deal with better

/**
 * A bloom filter that represents a window
 * within a larger bloom filter.
 */
class BloomFilterWindow : private BloomFilter
{
public:

	/** Constructor.
	 *
	 * @param fullBloomSize size in bits of the containing bloom filter
	 * @param startBitPos index of first bit in the window
	 * @param endBitPos index of last bit in the window
	 */
	BloomFilterWindow(size_t fullBloomSize, size_t startBitPos, size_t endBitPos) :
		BloomFilter(endBitPos - startBitPos + 1),
		m_fullBloomSize(fullBloomSize),
		m_startBitPos(startBitPos),
		m_endBitPos(endBitPos)
	{
		assert(startBitPos < fullBloomSize);
		assert(endBitPos < fullBloomSize);
		assert(startBitPos <= endBitPos);
	}

	/** Return the size of the bit array. */
	size_t size() const
	{
		return BloomFilter::size();
	}


	/** Return the number of elements with count >= MAX_COUNT. */
	size_t popcount() const
	{
		return BloomFilter::popcount();
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return BloomFilter::FPR();
	}

	/** Return whether the specified bit is set. */
	bool operator[](size_t i) const
	{
		if (i >= m_startBitPos && i <= m_endBitPos)
			return BloomFilter::operator[](i - m_startBitPos);
		return false;
	}

	/** Return whether the object is present in this set. */
	bool operator[](const Bloom::key_type& key) const
	{
		return (*this)[Bloom::hash(key) % m_fullBloomSize];
	}

	/** Add the object with the specified index to this set. */
	void insert(size_t i)
	{
		if (i >= m_startBitPos && i <= m_endBitPos)
			BloomFilter::insert(i - m_startBitPos);
	}

	/** Add the object to this set. */
	void insert(const Bloom::key_type& key)
	{
		insert(Bloom::hash(key) % m_fullBloomSize);
	}

	/**
	 * Return the full size of the containing bloom
	 * filter (in bits).
	 */
	size_t getFullBloomSize() const
	{
		return m_fullBloomSize;
	}

	/** Write the bloom filter to a stream */
	void write(std::ostream& out) const
	{
		out << BloomFilter::BLOOM_VERSION << '\n';
		assert(out);
		out << Kmer::length() << '\n';
		assert(out);
		writeBloomDimensions(out);
		assert(out);

		size_t bits = size();
		size_t bytes = (bits + 7) / 8;
		char buf[IO_BUFFER_SIZE];
		for (size_t i = 0, j = 0; i < bytes;) {
			size_t writeSize = std::min(IO_BUFFER_SIZE, bytes - i);
			for (size_t k = 0; k < writeSize; k++) {
				buf[k] = 0;
				for (unsigned l = 0; l < 8; l++, j++) {
					buf[k] <<= 1;
					if (j < bits && BloomFilter::m_array[j]) {
						buf[k] |= 1;
					}
				}
			}
			out.write(buf, writeSize);
			assert(out);
			i += writeSize;
		}
	}

	void read(std::istream& in, bool
			loadUnion = false, unsigned shrinkFactor = 1)
	{
		unsigned bloomVersion;
		in >> bloomVersion >> expect("\n");
		assert(in);
		if (bloomVersion != BLOOM_VERSION) {
			std::cerr << "error: bloom filter version (`"
				<< bloomVersion << "'), does not match version required "
				"by this program (`" << BLOOM_VERSION << "').\n";
			exit(EXIT_FAILURE);
		}

		// read bloom filter k value

		unsigned k;
		in >> k >> expect("\n");
		assert(in);
		if (k != Kmer::length()) {
			std::cerr << "error: this program must be run with the same kmer "
				"size as the bloom filter being loaded (k="
				<< k << ").\n";
			exit(EXIT_FAILURE);
		}

		// read bloom filter dimensions

		size_t size, startBitPos, endBitPos;
		readBloomDimensions(in, size, startBitPos, endBitPos);

		// shrink factor allows building a smaller
		// bloom filter from a larger one

		if (size % shrinkFactor != 0) {
			std::cerr << "error: the number of bits in the original bloom "
				"filter must be evenly divisible by the shrink factor (`"
				<< shrinkFactor << "')\n";
			exit(EXIT_FAILURE);
		}

		size /= shrinkFactor;

		if(loadUnion && size != this->size()) {
			std::cerr << "error: can't union bloom filters "
				"with different sizes.\n";
			exit(EXIT_FAILURE);
		} else {
			m_array.resize(size);
		}

		// read bit vector

		if (!loadUnion)
			m_array.reset();

		size_t offset = startBitPos;
		size_t bits = endBitPos - startBitPos + 1;
		size_t bytes = (bits + 7) / 8;

		char buf[IO_BUFFER_SIZE];
		for (size_t i = 0, j = offset; i < bytes; ) {
			size_t readSize = std::min(IO_BUFFER_SIZE, bytes - i);
			in.read(buf, readSize);
			assert(in);
			for (k = 0; k < readSize; k++) {
				for (unsigned l = 0; l < 8 && j < offset + bits; l++, j++) {
					bool bit = buf[k] & (1 << (7 - l));
					size_t index = j % size;
					m_array[index] |= bit;
				}
			}
			i += readSize;
		}
	}

	void loadSeq(unsigned k, const std::string& seq)
	{
		if (seq.size() < k)
			return;
		for (size_t i = 0; i < seq.size() - k + 1; ++i) {
			std::string kmer = seq.substr(i, k);
			size_t pos = kmer.find_last_not_of("ACGTacgt");
			if (pos == std::string::npos) {
				insert(Kmer(kmer));
			} else
				i += pos;
		}
	}

	/** Operator for writing the bloom filter to a stream */
	friend std::ostream& operator<<(std::ostream& out, const BloomFilterWindow& o)
	{
		o.write(out);
		return out;
	}

protected:

	void writeBloomDimensions(std::ostream& out) const
	{
		out << m_fullBloomSize
			<< '\t' << m_startBitPos
			<< '\t' << m_endBitPos
			<< '\n';
	}

	void readBloomDimensions(std::istream& in,
		size_t& size, size_t& startBitPos,
		size_t& endBitPos)
	{
		BloomFilter::readBloomDimensions(in, size,
				startBitPos, endBitPos);

		m_fullBloomSize = size;
		m_startBitPos = startBitPos;
		m_endBitPos = endBitPos;

		size = endBitPos - startBitPos + 1;
	}

private:

	size_t m_fullBloomSize;
	size_t m_startBitPos, m_endBitPos;
};

#endif
