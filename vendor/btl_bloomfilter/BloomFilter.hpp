/*
 *
 * BloomFilter.hpp
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_

#include "vendor/IOUtil.h"
#include "vendor/cpptoml/include/cpptoml.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace std;

static const uint8_t bitsPerChar = 0x08;
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

inline unsigned
popCnt(unsigned char x)
{
	return ((0x876543210 >> (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
	        ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2)) &
	       0xf;
}

class BloomFilter
{
  public:
	/*
	 * Default constructor.
	 */
	BloomFilter()
	  : m_filter(NULL)
	  , m_size(0)
	  , m_sizeInBytes(0)
	  , m_hashNum(0)
	  , m_kmerSize(0)
	  , m_dFPR(0)
	  , m_nEntry(0)
	  , m_tEntry(0)
	{}

	/* De novo filter constructor.
	 *
	 * preconditions:
	 * filterSize must be a multiple of 64
	 *
	 * kmerSize refers to the number of bases the kmer has
	 */
	BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize)
	  : m_filter(NULL)
	  , m_size(filterSize)
	  , m_hashNum(hashNum)
	  , m_kmerSize(kmerSize)
	  , m_dFPR(0)
	  , m_nEntry(0)
	  , m_tEntry(0)
	{
		initSize(m_size);
	}

	/* De novo filter constructor.
	 * Allocates a filter size based on the number of expected elements and FPR
	 *
	 * If hashNum is set to 0, an optimal value is computed based on the FPR
	 */
	BloomFilter(size_t expectedElemNum, double fpr, unsigned hashNum, unsigned kmerSize)
	  : m_size(0)
	  , m_hashNum(hashNum)
	  , m_kmerSize(kmerSize)
	  , m_dFPR(fpr)
	  , m_nEntry(0)
	  , m_tEntry(0)
	{
		if (m_hashNum == 0) {
			m_hashNum = calcOptiHashNum(m_dFPR);
		}
		if (m_size == 0) {
			m_size = calcOptimalSize(expectedElemNum, m_dFPR);
		}
		initSize(m_size);
	}

	BloomFilter(const string& filterFilePath)
	  : m_filter(NULL)
	{
		loadFilter(filterFilePath);
	}

	void loadFilter(const string& filterFilePath)
	{
		std::ifstream file(filterFilePath);
		assert_good(file, filterFilePath);
		loadHeader(file);
		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
		file.read(reinterpret_cast<char*>(m_filter), m_sizeInBytes);
		assert_good(file, filterFilePath);
		file.close();
	}

	void loadHeader(std::istream& file)
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
		m_dFPR = *bloomFilterTable->get_as<double>("dFPR");
		m_nEntry = *bloomFilterTable->get_as<uint64_t>("nEntry");
		m_tEntry = *bloomFilterTable->get_as<uint64_t>("Entry");
		initSize(m_size);
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	void insert(vector<uint64_t> const& precomputed)
	{

		// iterates through hashed values adding it to the filter
		for (unsigned i = 0; i < m_hashNum; ++i) {
			uint64_t normalizedValue = precomputed.at(i) % m_size;
			__sync_or_and_fetch(
			    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	void insert(const uint64_t precomputed[])
	{

		// iterates through hashed values adding it to the filter
		for (unsigned i = 0; i < m_hashNum; ++i) {
			uint64_t normalizedValue = precomputed[i] % m_size;
			__sync_or_and_fetch(
			    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * Returns if already inserted
	 */
	bool insertAndCheck(const uint64_t precomputed[])
	{
		// iterates through hashed values adding it to the filter
		bool found = true;
		for (unsigned i = 0; i < m_hashNum; ++i) {
			uint64_t normalizedValue = precomputed[i] % m_size;
			found &= __sync_fetch_and_or(
			             &m_filter[normalizedValue / bitsPerChar],
			             bitMask[normalizedValue % bitsPerChar]) >>
			             (normalizedValue % bitsPerChar) &
			         1;
		}
		return found;
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 * Returns if already inserted
	 */
	bool insertAndCheck(vector<uint64_t> const& precomputed)
	{
		// iterates through hashed values adding it to the filter
		bool found = true;
		for (unsigned i = 0; i < m_hashNum; ++i) {
			uint64_t normalizedValue = precomputed.at(i) % m_size;
			found &= __sync_fetch_and_or(
			             &m_filter[normalizedValue / bitsPerChar],
			             bitMask[normalizedValue % bitsPerChar]) >>
			             (normalizedValue % bitsPerChar) &
			         1;
		}
		return found;
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	bool contains(vector<uint64_t> const& precomputed) const
	{
		for (unsigned i = 0; i < m_hashNum; ++i) {
			uint64_t normalizedValue = precomputed.at(i) % m_size;
			unsigned char bit = bitMask[normalizedValue % bitsPerChar];
			if ((m_filter[normalizedValue / bitsPerChar] & bit) != bit) {
				return false;
			}
		}
		return true;
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	bool contains(const uint64_t precomputed[]) const
	{
		for (unsigned i = 0; i < m_hashNum; ++i) {
			uint64_t normalizedValue = precomputed[i] % m_size;
			unsigned char bit = bitMask[normalizedValue % bitsPerChar];
			if ((m_filter[normalizedValue / bitsPerChar] & bit) != bit) {
				return false;
			}
		}
		return true;
	}

	void writeHeader(std::ostream& out) const
	{
		/* Initialize cpptoml root table
		   Note: Tables and fields are unordered
		   Ordering of table is maintained by directing the table
		   to the output stream immediately after completion  */
		std::shared_ptr<cpptoml::table> root = cpptoml::make_table();

		/* Initialize bloom filter section and insert fields
		   and output to ostream */
		auto header = cpptoml::make_table();
		header->insert("KmerSize", m_kmerSize);
		header->insert("HashNum", m_hashNum);
		header->insert("BloomFilterSize", m_size);
		header->insert("BloomFilterSizeInBytes", m_sizeInBytes);
		header->insert("dFPR", m_dFPR);
		header->insert("nEntry", m_nEntry);
		header->insert("Entry", m_tEntry);
		std::string magic(MAGIC_HEADER_STRING);
		root->insert(magic, header);
		out << *root;

		// Output [HeaderEnd]\n to ostream to mark the end of the header
		out << "[HeaderEnd]\n";
	}

	/** Serialize the Bloom filter to a stream */
	friend std::ostream& operator<<(std::ostream& out, const BloomFilter& bloom)
	{
		bloom.writeHeader(out);
		// NOLINTNEXTLINE(google-readability-casting)
		out.write(reinterpret_cast<char*>(bloom.m_filter), bloom.m_sizeInBytes);
		return out;
	}

	/*
	 * Stores the filter as a binary file to the path specified
	 * Stores uncompressed because the random data tends to
	 * compress poorly anyway
	 */
	void storeFilter(const string& filterFilePath) const
	{
		std::ofstream ofs(filterFilePath.c_str(), std::ios::out | std::ios::binary);
		assert_good(ofs, filterFilePath);
		std::cerr << "Writing a " << m_sizeInBytes << " byte filter to " << filterFilePath
		          << " on disk.\n";
		ofs << *this;
		ofs.flush();
		assert_good(ofs, filterFilePath);
		ofs.close();
	}

	uint64_t getPop() const
	{
		uint64_t i, popBF = 0;
		//#pragma omp parallel for reduction(+:popBF)
		for (i = 0; i < (m_size + 7) / 8; i++)
			popBF = popBF + popCnt(m_filter[i]);
		return popBF;
	}

	unsigned getHashNum() const { return m_hashNum; }

	unsigned getKmerSize() const { return m_kmerSize; }

	/*
	 * Calculates that False positive rate that a redundant entry is actually
	 * a unique entry
	 */
	double getRedudancyFPR()
	{
		assert(m_nEntry > 0);
		double total = log(calcFPR_numInserted(1));
		for (uint64_t i = 2; i < m_nEntry; ++i) {
			total = log(exp(total) + calcFPR_numInserted(i));
		}
		return exp(total) / m_nEntry;
	}

	/*
	 * Return FPR based on popcount
	 */
	double getFPR() const { return pow(double(getPop()) / double(m_size), double(m_hashNum)); }

	/*
	 * Return FPR based on number of inserted elements
	 */
	double getFPR_numEle() const
	{
		assert(m_nEntry > 0);
		return calcFPR_numInserted(m_nEntry);
	}

	uint64_t getnEntry() { return m_nEntry; }

	uint64_t gettEntry() { return m_tEntry; }

	void setnEntry(uint64_t value) { m_nEntry = value; }

	void settEntry(uint64_t value) { m_tEntry = value; }

	uint64_t getFilterSize() const { return m_size; }

	uint64_t sizeInBytes() const { return m_sizeInBytes; }

	~BloomFilter() { delete[] m_filter; }

  protected:
	BloomFilter(const BloomFilter& that); // to prevent copy construction

	/*
	 * Checks filter size and initializes filter
	 */
	void initSize(size_t size)
	{
		if (size % 8 != 0) {
			cerr << "ERROR: Filter Size \"" << size << "\" is not a multiple of 8." << endl;
			exit(1);
		}
		m_sizeInBytes = size / bitsPerChar;
		if (m_filter != NULL)
			delete[] m_filter;
		m_filter = new unsigned char[m_sizeInBytes]();
	}

	/*
	 * Only returns multiples of 64 for filter building purposes
	 * Is an estimated size using approximations of FPR formula
	 * given the number of hash functions
	 */
	size_t calcOptimalSize(size_t entries, double fpr) const
	{
		size_t non64ApproxVal = size_t(
		    -double(entries) * double(m_hashNum) /
		    log(1.0 - pow(fpr, double(1 / double(m_hashNum)))));

		return non64ApproxVal + (64 - non64ApproxVal % 64);
	}

	/*
	 * Calculates the optimal number of hash function to use
	 * Calculation assumes optimal ratio of bytes per entry given a fpr
	 */
	static unsigned calcOptiHashNum(double fpr) { return unsigned(-log(fpr) / log(2)); }

	/*
	 * Calculate FPR based on hash functions, size and number of entries
	 * see http://en.wikipedia.org/wiki/Bloom_filter
	 */
	double calcFPR_numInserted(size_t numEntr) const
	{
		return pow(
		    1.0 - pow(1.0 - 1.0 / double(m_size), double(numEntr) * m_hashNum), double(m_hashNum));
	}

	/*
	 * Calculates the optimal FPR to use based on hash functions
	 */
	double calcFPR_hashNum(unsigned hashFunctNum) const { return pow(2, -double(hashFunctNum)); }

	uint8_t* m_filter;
	size_t m_size;
	size_t m_sizeInBytes;
	unsigned m_hashNum;
	unsigned m_kmerSize;
	double m_dFPR;
	uint64_t m_nEntry;
	uint64_t m_tEntry;
	static constexpr const char* MAGIC_HEADER_STRING = "BTLBloomFilter_v1";
};

#endif /* BLOOMFILTER_H_ */
