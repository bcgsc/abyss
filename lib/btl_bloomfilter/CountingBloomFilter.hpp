// Written by Johnathan Wong and Sauparna Palchowdhury.

#ifndef COUNTINGBLOOMFILTER_HPP // NOLINT(llvm-header-guard)
#define COUNTINGBLOOMFILTER_HPP

#include "vendor/IOUtil.h"
#include "vendor/cpptoml/include/cpptoml.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

// Forward declaraions.
template<typename T>
class CountingBloomFilter;

// Method declarations.
template<typename T>
std::ostream&
operator<<(std::ostream& out, const CountingBloomFilter<T>& bloom);

template<typename T>
class CountingBloomFilter
{
  public:
	CountingBloomFilter() = default;
	CountingBloomFilter(
	    size_t sizeInBytes,
	    unsigned hashNum,
	    unsigned kmerSize,
	    unsigned countThreshold)
	  : m_hashNum(hashNum)
	  , m_kmerSize(kmerSize)
	  , m_countThreshold(countThreshold)
	{
		int remainder = sizeInBytes % 8;
		if (remainder == 0) {
			m_sizeInBytes = sizeInBytes;
			m_size = m_sizeInBytes / sizeof(T);
			m_filter.resize(m_size, 0);
		} else {
			m_sizeInBytes = sizeInBytes + 8 - remainder;
			m_size = m_sizeInBytes / sizeof(T);
			m_filter.resize(m_size, 0);
		}
	}
	CountingBloomFilter(const std::string& path, unsigned countThreshold);
	T operator[](size_t i) { return m_filter[i]; }
	template<typename U>
	T minCount(const U& hashes) const
	{
		T min = m_filter[hashes[0] % m_size];
		for (size_t i = 1; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_size;
			if (m_filter[pos] < min) {
				min = m_filter[pos];
			}
		}
		return min;
	}
	template<typename U>
	bool contains(const U& hashes) const;
	template<typename U>
	void insert(const U& hashes);
	template<typename U>
	bool insertAndCheck(const U& hashes);
	template<typename U>
	void incrementMin(const U& hashes);
	template<typename U>
	void incrementAll(const U& hashes);
	unsigned getKmerSize() const { return m_kmerSize; };
	unsigned getHashNum() const { return m_hashNum; };
	unsigned threshold() const { return m_countThreshold; };
	size_t size() const { return m_size; };
	size_t sizeInBytes() const { return m_sizeInBytes; };
	size_t popCount() const;
	size_t filtered_popcount() const;
	double FPR() const;
	double filtered_FPR() const;
	void loadHeader(std::istream& file);
	void loadFilter(const std::string& path);
	void storeHeader(std::ostream& out) const;
	void storeFilter(const std::string& path) const;
	friend std::ostream& operator<<<>(std::ostream&, const CountingBloomFilter&);

  private:
	// m_filter             : A vector of elements of type T.
	// m_size               : Size of bloom filter (size of m_filter array).
	// m_sizeInBytes        : Size of the bloom filter in bytes, that is,
	//                        (m_size * sizeof(T)).
	// m_hashNum            : Number of hash functions.
	// m_kmerSize           : Size of a k-mer.
	// m_countThreshold     : A count greater or equal to this threshold
	//                        establishes existence of an element in the filter.
	// m_bitsPerCounter     : Number of bits per counter.
	// MAGIC_HEADER_STRING  : Magic string used to identify the type of bloom filter.

	std::vector<T> m_filter;
	size_t m_size = 0;
	size_t m_sizeInBytes = 0;
	unsigned m_hashNum = 0;
	unsigned m_kmerSize = 0;
	unsigned m_countThreshold = 0;
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers, readability-magic-numbers)
	unsigned m_bitsPerCounter = 8;
	static constexpr const char* MAGIC_HEADER_STRING = "BTLCountingBloomFilter_v1";
};

// Method definitions

// Bloom filter operations

/*
  Use of atomic increments in incrementMin() and incrementAll():

  A atomic compare-and-swap (CAS) operation increments m_filter[pos]. The CAS
  operation takes a memory location and a value that the caller believes the
  location currently stores. If the memory location still holds that value when
  the atomic compare-and-swap executes, then a new value is stored and 'true'
  returned; otherwise, memory is left unchanged and 'false' returned.

  The value of m_filter[pos] may be changed by another thread between a read from
  that memory location and a write to it. The CAS operation is called in a loop
  until it succeeds, which ensures that a write does not happen if some other
  thread has incremented the value between this thread's read and write.

  Note that CAS operations suffer from the ABA problem.
*/

// Of the m_hashNum counters, increment all the minimum values.
template<typename T>
template<typename U>
inline void
CountingBloomFilter<T>::incrementMin(const U& hashes)
{
	// update flag to track if increment is done on at least one counter
	bool updateDone = false;
	T newVal;
	T minVal = minCount(hashes);
	while (!updateDone) {
		// Simple check to deal with overflow
		newVal = minVal + 1;
		if (minVal > newVal) {
			return;
		}
		for (size_t i = 0; i < m_hashNum; ++i) {
			// NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg, hicpp-vararg)
			if (__sync_bool_compare_and_swap(&m_filter[hashes[i] % m_size], minVal, newVal)) {
				updateDone = true;
			}
		}
		// Recalculate minval because if increment fails, it needs a new minval to use and
		// if it doesnt hava a new one, the while loop runs forever.
		if (!updateDone) {
			minVal = minCount(hashes);
		}
	}
}

// Increment all the m_hashNum counters.
template<typename T>
template<typename U>
inline void
CountingBloomFilter<T>::incrementAll(const U& hashes)
{
	T currentVal;
	T newVal;
	for (size_t i = 0; i < m_hashNum; ++i) {
		size_t pos = hashes[i] % m_size;
		do {
			currentVal = m_filter[pos];
			newVal = currentVal + 1;
			if (newVal < currentVal) {
				break;
			}
			// NOLINTNEXTLINE(cppcoreguidelines-pro-type-vararg, hicpp-vararg)
		} while (!__sync_bool_compare_and_swap(&m_filter[pos], currentVal, newVal));
	}
}

// Check if an element exists. If the minimum count at the m_hashNum positions
// of m_filter is more than or equal to a predefined count threshold, then the
// element is said to be present in the Bloom filter. count() therefore returns
// true when this condition is satisfied, or else, false.

template<typename T>
template<typename U>
inline bool
CountingBloomFilter<T>::contains(const U& hashes) const
{
	return minCount(hashes) >= m_countThreshold;
}

template<typename T>
template<typename U>
inline void
CountingBloomFilter<T>::insert(const U& hashes)
{
	incrementMin(hashes);
}

template<typename T>
template<typename U>
inline bool
CountingBloomFilter<T>::insertAndCheck(const U& hashes)
{
	bool found = contains(hashes);
	incrementMin(hashes);
	return found;
}

/* Count the number of non-zero counters. */
template<typename T>
size_t
CountingBloomFilter<T>::popCount() const
{
	size_t count = 0;
	for (size_t i = 0; i < m_size; ++i) {
		if (m_filter[i] != 0) {
			++count;
		}
	}
	return count;
}

/* Count the number of above threshold counters. */
template<typename T>
size_t
CountingBloomFilter<T>::filtered_popcount() const
{
	size_t count = 0;
	for (size_t i = 0; i < m_size; ++i) {
		if (m_filter[i] >= m_countThreshold) {
			++count;
		}
	}
	return count;
}

template<typename T>
double
CountingBloomFilter<T>::FPR() const
{
	// NOLINTNEXTLINE(google-readability-casting)
	return std::pow((double)popCount() / (double)m_size, m_hashNum);
}

template<typename T>
double
CountingBloomFilter<T>::filtered_FPR() const
{
	// NOLINTNEXTLINE(google-readability-casting)
	return std::pow((double)filtered_popcount() / (double)m_size, m_hashNum);
}

// Serialization interface.
template<typename T>
CountingBloomFilter<T>::CountingBloomFilter(const std::string& path, unsigned countThreshold)
  : m_countThreshold(countThreshold)
{
	loadFilter(path);
}

template<typename T>
void
CountingBloomFilter<T>::loadFilter(const std::string& path)
{
	std::ifstream file(path);
	assert_good(file, path);
	loadHeader(file);
	m_filter.resize(m_sizeInBytes);
	// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
	file.read(reinterpret_cast<char*>(m_filter.data()), m_sizeInBytes);
	assert_good(file, path);
	file.close();
}

template<typename T>
void
CountingBloomFilter<T>::loadHeader(std::istream& file)
{
	std::string magic_header(MAGIC_HEADER_STRING);
	(magic_header.insert(0, "[")).append("]");
	std::string line;
	std::getline(file, line);
	if (line != magic_header) {
		std::cerr << "ERROR: magic string does not match (likely version mismatch)\n"
		          << "Your magic string:                " << line << "\n"
		          << "CountingBloomFilter magic string: " << magic_header << std::endl;
		exit(EXIT_FAILURE);
	}

	/* Read bloom filter line by line until it sees "[HeaderEnd]"
	   which is used to mark the end of the header section and
	   assigns the header to a char array*/
	std::string headerEnd = "[HeaderEnd]";
	std::string toml_buffer((line + "\n"));
	bool headerEndCheck = false;
	while (std::getline(file, line)) {
		toml_buffer.append(line + "\n");
		if (line == headerEnd) {
			headerEndCheck = true;
			break;
		}
	}
	if (!headerEndCheck) {
		std::cerr << "ERROR: pre-built bloom filter does not have the correct header end."
		          << std::endl;
		exit(EXIT_FAILURE);
	}

	// Send the char array to a stringstream for the cpptoml parser to parse
	std::istringstream toml_stream(toml_buffer);
	cpptoml::parser toml_parser(toml_stream);
	auto header_config = toml_parser.parse();

	// Obtain header values from toml parser and assign them to class members
	std::string magic(MAGIC_HEADER_STRING);
	auto bloomFilterTable = header_config->get_table(magic);
	m_size = *bloomFilterTable->get_as<size_t>("BloomFilterSize");
	m_hashNum = *bloomFilterTable->get_as<unsigned>("HashNum");
	m_kmerSize = *bloomFilterTable->get_as<unsigned>("KmerSize");
	m_sizeInBytes = *bloomFilterTable->get_as<size_t>("BloomFilterSizeInBytes");
	m_bitsPerCounter = *bloomFilterTable->get_as<unsigned>("BitsPerCounter");
}

template<typename T>
void
CountingBloomFilter<T>::storeFilter(const std::string& path) const
{
	std::ofstream ofs(path.c_str(), std::ios::out | std::ios::binary);
	assert_good(ofs, path);
	std::cerr << "Writing a " << m_sizeInBytes << " byte filter to " << path << " on disk.\n";
	ofs << *this;
	ofs.flush();
	assert_good(ofs, path);
	ofs.close();
}

template<typename T>
void
CountingBloomFilter<T>::storeHeader(std::ostream& out) const
{
	/* Initialize cpptoml root table
	   Note: Tables and fields are unordered
	   Ordering of table is maintained by directing the table
	   to the output stream immediately after completion  */
	std::shared_ptr<cpptoml::table> root = cpptoml::make_table();

	/* Initialize bloom filter section and insert fields
	   and output to ostream */
	auto header = cpptoml::make_table();
	header->insert("BitsPerCounter", m_bitsPerCounter);
	header->insert("KmerSize", m_kmerSize);
	header->insert("HashNum", m_hashNum);
	header->insert("BloomFilterSize", m_size);
	header->insert("BloomFilterSizeInBytes", m_sizeInBytes);
	std::string magic(MAGIC_HEADER_STRING);
	root->insert(magic, header);
	out << *root;

	// Output [HeaderEnd]\n to ostream to mark the end of the header
	out << "[HeaderEnd]\n";
}

// Serialize the bloom filter to a C++ stream
template<typename T>
std::ostream&
operator<<(std::ostream& out, const CountingBloomFilter<T>& bloom)
{
	bloom.storeHeader(out);
	// NOLINTNEXTLINE(google-readability-casting)
	out.write((const char*)bloom.m_filter.data(), bloom.m_sizeInBytes);
	return out;
}

#endif // COUNTINGBLOOMFILTER_HPP
