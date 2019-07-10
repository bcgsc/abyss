#ifndef COUNTINGBLOOM_H
#define COUNTINGBLOOM_H

#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "cpptoml/include/cpptoml.h"
#include "IOUtil.h"

// Forward declaraions.
template<typename T>
class CountingBloomFilter;

// Method declarations.
template<typename T>
std::ostream&
operator<<(std::ostream& os, const CountingBloomFilter<T>& cbf);

template<typename T>
class CountingBloomFilter
{
  public:
	CountingBloomFilter() {}
	CountingBloomFilter(size_t sz, unsigned hashNum, unsigned kmerSize, unsigned countThreshold)
	  : m_filter(new T[sz])
	  , m_size(sz)
	  , m_sizeInBytes(sz * sizeof(T))
	  , m_hashNum(hashNum)
	  , m_kmerSize(kmerSize)
	  , m_countThreshold(countThreshold)
	{
		std::memset(m_filter, 0, m_sizeInBytes);
	}
	CountingBloomFilter(const std::string& path, unsigned countThreshold);
	~CountingBloomFilter() { delete[] m_filter; }
	T operator[](size_t i) { return m_filter[i]; }
	template<typename U>
	T minCount(const U& hashes) const
	{
		T min = m_filter[hashes[0] % m_size];
		for (size_t i = 1; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_size;
			if (m_filter[pos] < min)
				min = m_filter[pos];
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
	void readHeader(std::istream& file);
	void readFilter(const std::string& path);
	void writeHeader(std::ostream& out) const;
	void writeFilter(const std::string& path) const;
	friend std::ostream& operator<<<>(std::ostream&, const CountingBloomFilter&);

  private:
	// m_filter             : An array of elements of type T; the bit-array or
	//                        filter.
	// m_size               : Size of bloom filter (size of m_filter array).
	// m_sizeInBytes        : Size of the bloom filter in bytes, that is,
	//                        (m_size * sizeof(T)).
	// m_hashNum            : Number of hash functions.
	// m_kmerSize           : Size of a k-mer.
	// m_countThreshold     : A count greater or equal to this threshold
	//                        establishes existence of an element in the filter.
	// m_bitsPerCounter     : Number of bits per counter.
	// MAGIC_HEADER_STRING  : Magic string used to identify the type of bloom filter.

	T* m_filter = nullptr;
	size_t m_size = 0;
	size_t m_sizeInBytes = 0;
	unsigned m_hashNum = 0;
	unsigned m_kmerSize = 0;
	unsigned m_countThreshold = 0;
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
	return;
}

// Increment all the m_hashNum counters.
template<typename T>
template<typename U>
inline void
CountingBloomFilter<T>::incrementAll(const U& hashes)
{
	T currentVal, newVal;
	for (size_t i = 0; i < m_hashNum; ++i) {
		size_t pos = hashes[i] % m_size;
		do {
			currentVal = m_filter[pos];
			newVal = currentVal + 1;
			if (newVal < currentVal) {
				break;
			}
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
	return std::pow((double)popCount() / (double)m_size, m_hashNum);
}

template<typename T>
double
CountingBloomFilter<T>::filtered_FPR() const
{
	return std::pow((double)filtered_popcount() / (double)m_size, m_hashNum);
}

// Serialization interface.
template<typename T>
CountingBloomFilter<T>::CountingBloomFilter(const std::string& path, unsigned countThreshold)
  : m_countThreshold(countThreshold)
{
	readFilter(path);
}

template<typename T>
void
CountingBloomFilter<T>::readFilter(const std::string& path)
{
	std::ifstream file(path);
	if (!file) {
		std::cerr << "ERROR: Failed to open file: " << path << "\n";
		exit(EXIT_FAILURE);
	}
	readHeader(file);
	char* filter = new char[m_sizeInBytes]; 
	file.read(filter, m_sizeInBytes);
	m_filter = reinterpret_cast<T*>(filter);
	if (!file) {
		std::cerr << "ERROR: The byte array could not be read from the file: " << path << "\n";
		exit(EXIT_FAILURE);
	}
	file.close();
}

template<typename T>
void
CountingBloomFilter<T>::readHeader(std::istream& file)
{
	std::string magic_header(MAGIC_HEADER_STRING);
	(magic_header.insert(0, "[")).append("]");
	std::string line;
	std::getline(file, line);
	if (line.compare(magic_header) != 0) {
		std::cerr << "ERROR: magic string does not match (likely version mismatch)\n"
		          << "Your magic string:                " << line << "\n"
		          << "CountingBloomFilter magic string: " << magic_header << std::endl;
		exit(EXIT_FAILURE);
	}

	/* Read bloom filter line by line until it sees "[HeaderEnd]"
	   which is used to mark the end of the header section and
	   assigns the header to a char array*/
	std::string headerEnd = "[HeaderEnd]";
	char* toml_buffer = new char[0];
	while (std::getline(file, line)) {
		if (line == headerEnd) {
			int currPos = file.tellg();
			delete[] toml_buffer;
			toml_buffer = new char[currPos];
			file.seekg(0, file.beg);
			file.read(toml_buffer, currPos);
			file.seekg(currPos, file.beg);
			break;
		}
	}

	// Send the char array to a stringstream for the cpptoml parser to parse
	std::istringstream toml_stream(toml_buffer);
	delete[] toml_buffer;
	cpptoml::parser toml_parser{ toml_stream };
	auto header_config = toml_parser.parse();

	// Obtain header values from toml parser and assign them to class members
	std::string magic(MAGIC_HEADER_STRING);
	auto bloomFilterTable = header_config->get_table(magic);
	auto toml_size = bloomFilterTable->get_as<size_t>("BloomFilterSize");
	auto toml_kmerSize = bloomFilterTable->get_as<unsigned>("KmerSize");
	auto toml_bitsPerCounter = bloomFilterTable->get_as<unsigned>("BitsPerCounter");
	auto toml_hashNum = bloomFilterTable->get_as<unsigned>("HashNum");
	m_size = *toml_size;
	m_hashNum = *toml_hashNum;
	m_kmerSize = *toml_kmerSize;
	m_sizeInBytes = m_size * sizeof(T);
	m_bitsPerCounter = *toml_bitsPerCounter;
}

template<typename T>
void
CountingBloomFilter<T>::writeFilter(const std::string& path) const
{
	std::ofstream ofs(path.c_str(), std::ios::out | std::ios::binary);
	std::cerr << "Writing a " << m_sizeInBytes << " byte filter to a file on disk.\n";
	ofs << *this;
	ofs.close();
	assert_good(ofs, path);
}

template<typename T>
void
CountingBloomFilter<T>::writeHeader(std::ostream& out) const
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
	std::string magic(MAGIC_HEADER_STRING);
	root->insert(magic, header);
	out << *root;

	/* Initalize new cpptoml root table and HeaderEnd section,
	   and output to ostream */
	root = cpptoml::make_table();
	auto ender = cpptoml::make_table();
	root->insert(std::string("HeaderEnd"), ender);
	out << (*root);
	assert(out);
}

// Serialize the bloom filter to a C++ stream
template<typename T>
std::ostream&
operator<<(std::ostream& out, const CountingBloomFilter<T>& bloom)
{
	assert(out);
	bloom.writeHeader(out);
	assert(out);
	out.write(reinterpret_cast<char*>(bloom.m_filter), bloom.m_sizeInBytes);
	assert(out);
	return out;
}

#endif // COUNTINGBLOOM_H
