/*
 *
 * BloomFilter_pythonwrapper.cpp
 *
 *  Created on: Dec 28, 2015
 *
 */

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_
#include "rolling.h"
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
	{}

	/* De novo filter constructor.
	 *
	 * preconditions:
	 * filterSize must be a multiple of 64
	 * kmerSize refers to the number of bases the kmer has
	 * k-mers supplied to this object should be binary (2 bits per base)
	 */
	BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize)
	  : m_size(filterSize)
	  , m_hashNum(hashNum)
	  , m_kmerSize(kmerSize)
	{
		initSize(m_size);
		memset(m_filter, 0, m_sizeInBytes);
	}

	/*
	 * Loads the filter (file is a .bf file) from path specified
	 */
	BloomFilter(
	    size_t filterSize,
	    unsigned hashNum,
	    unsigned kmerSize,
	    string const& filterFilePath)
	  : m_size(filterSize)
	  , m_hashNum(hashNum)
	  , m_kmerSize(kmerSize)
	{
		initSize(m_size);

		FILE* file = fopen(filterFilePath.c_str(), "rb");
		if (file == NULL) {
			cerr << "file \"" << filterFilePath << "\" could not be read." << endl;
			exit(1);
		}

		long int lCurPos = ftell(file);
		fseek(file, 0, 2);
		size_t fileSize = ftell(file);
		fseek(file, lCurPos, 0);
		if (fileSize != m_sizeInBytes) {
			cerr << "Error: " << filterFilePath
			     << " does not match size given by its information file. Size: " << fileSize
			     << " vs " << m_sizeInBytes << " bytes." << endl;
			exit(1);
		}

		size_t countRead = fread(m_filter, fileSize, 1, file);
		if (countRead != 1 && fclose(file) != 0) {
			cerr << "file \"" << filterFilePath << "\" could not be read." << endl;
			exit(1);
		}
	}

	/*
	 * For precomputing hash values. kmerSize is the number of bytes of the original string used.
	 */
	vector<size_t> multiHash(const char* kmer) const
	{
		vector<size_t> tempHashValues(m_hashNum);
		uint64_t hVal = getChval(kmer, m_kmerSize);
		for (size_t i = 0; i < m_hashNum; ++i) {
			tempHashValues[i] = (rol(varSeed, i) ^ hVal);
		}
		return tempHashValues;
	}

	/*
	 * For precomputing hash values. kmerSize is the number of bytes of the original string used.
	 */
	vector<size_t> multiHash(const char* kmer, uint64_t& fhVal, uint64_t& rhVal) const
	{
		vector<size_t> tempHashValues(m_hashNum);
		fhVal = getFhval(kmer, m_kmerSize);
		rhVal = getRhval(kmer, m_kmerSize);
		uint64_t hVal = (rhVal < fhVal) ? rhVal : fhVal;
		for (unsigned i = 0; i < m_hashNum; i++) {
			tempHashValues[i] = (rol(varSeed, i) ^ hVal);
		}
		return tempHashValues;
	}

	/*
	 * For precomputing hash values. kmerSize is the number of bytes of the original string used.
	 */
	vector<size_t>
	multiHash(uint64_t& fhVal, uint64_t& rhVal, const char charOut, const char charIn) const
	{
		vector<size_t> tempHashValues(m_hashNum);
		fhVal = rol(fhVal, 1) ^ rol(seedTab[(unsigned char)charOut], m_kmerSize) ^
		        seedTab[(unsigned char)charIn];
		rhVal = ror(rhVal, 1) ^ ror(seedTab[(unsigned char)(charOut + cpOff)], 1) ^
		        rol(seedTab[(unsigned char)(charIn + cpOff)], m_kmerSize - 1);
		uint64_t hVal = (rhVal < fhVal) ? rhVal : fhVal;
		for (unsigned i = 0; i < m_hashNum; i++) {
			tempHashValues[i] = (rol(varSeed, i) ^ hVal);
		}
		return tempHashValues;
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	void insert(vector<size_t> const& precomputed)
	{

		// iterates through hashed values adding it to the filter
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t normalizedValue = precomputed.at(i) % m_size;
			__sync_or_and_fetch(
			    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
			//		m_filter[normalizedValue / bitsPerChar] |= bitMask[normalizedValue
			//				% bitsPerChar];
		}
	}

	void insert(const char* kmer)
	{
		uint64_t hVal = getChval(kmer, m_kmerSize);
		for (unsigned i = 0; i < m_hashNum; i++) {
			size_t normalizedValue = (rol(varSeed, i) ^ hVal) % m_size;
			__sync_or_and_fetch(
			    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
		}
	}

	void insert(const char* kmer, uint64_t& fhVal, uint64_t& rhVal)
	{
		fhVal = getFhval(kmer, m_kmerSize);
		rhVal = getRhval(kmer, m_kmerSize);
		uint64_t hVal = (rhVal < fhVal) ? rhVal : fhVal;
		for (unsigned i = 0; i < m_hashNum; i++) {
			size_t normalizedValue = (rol(varSeed, i) ^ hVal) % m_size;
			__sync_or_and_fetch(
			    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
		}
	}

	void insert(uint64_t& fhVal, uint64_t& rhVal, const char charOut, const char charIn)
	{
		fhVal = rol(fhVal, 1) ^ rol(seedTab[(unsigned char)charOut], m_kmerSize) ^
		        seedTab[(unsigned char)charIn];
		rhVal = ror(rhVal, 1) ^ ror(seedTab[(unsigned char)(charOut + cpOff)], 1) ^
		        rol(seedTab[(unsigned char)(charIn + cpOff)], m_kmerSize - 1);
		uint64_t hVal = (rhVal < fhVal) ? rhVal : fhVal;
		for (unsigned i = 0; i < m_hashNum; i++) {
			size_t normalizedValue = (rol(varSeed, i) ^ hVal) % m_size;
			__sync_or_and_fetch(
			    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	bool contains(vector<size_t> const& values)
	{
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t normalizedValue = values.at(i) % m_size;
			unsigned char bit = bitMask[normalizedValue % bitsPerChar];
			if ((m_filter[normalizedValue / bitsPerChar] & bit) != bit) {
				return false;
			}
		}
		return true;
	}

	/*
	 * Single pass filtering, computes hash values on the fly
	 */
	bool contains(const char* kmer)
	{
		uint64_t hVal = getChval(kmer, m_kmerSize);
		for (unsigned i = 0; i < m_hashNum; i++) {
			size_t normalizedValue = (rol(varSeed, i) ^ hVal) % m_size;
			unsigned char bit = bitMask[normalizedValue % bitsPerChar];
			if ((m_filter[normalizedValue / bitsPerChar] & bit) == 0)
				return false;
		}
		return true;
	}

	bool contains(const char* kmer, uint64_t& fhVal, uint64_t& rhVal)
	{
		fhVal = getFhval(kmer, m_kmerSize);
		rhVal = getRhval(kmer, m_kmerSize);
		uint64_t hVal = (rhVal < fhVal) ? rhVal : fhVal;
		for (unsigned i = 0; i < m_hashNum; i++) {
			size_t normalizedValue = (rol(varSeed, i) ^ hVal) % m_size;
			unsigned char bit = bitMask[normalizedValue % bitsPerChar];
			if ((m_filter[normalizedValue / bitsPerChar] & bit) == 0)
				return false;
		}
		return true;
	}

	bool contains(uint64_t& fhVal, uint64_t& rhVal, const char charOut, const char charIn)
	{
		fhVal = rol(fhVal, 1) ^ rol(seedTab[(unsigned char)charOut], m_kmerSize) ^
		        seedTab[(unsigned char)charIn];
		rhVal = ror(rhVal, 1) ^ ror(seedTab[(unsigned char)(charOut + cpOff)], 1) ^
		        rol(seedTab[(unsigned char)(charIn + cpOff)], m_kmerSize - 1);
		uint64_t hVal = (rhVal < fhVal) ? rhVal : fhVal;
		for (unsigned i = 0; i < m_hashNum; i++) {
			size_t normalizedValue = (rol(varSeed, i) ^ hVal) % m_size;
			unsigned char bit = bitMask[normalizedValue % bitsPerChar];
			if ((m_filter[normalizedValue / bitsPerChar] & bit) == 0)
				return false;
		}
		return true;
	}

	/*
	 * Stores the filter as a binary file to the path specified
	 * Stores uncompressed because the random data tends to
	 * compress poorly anyway
	 */
	void storeFilter(string const& filterFilePath) const
	{
		ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

		cerr << "Storing filter. Filter is " << m_sizeInBytes << "bytes." << endl;

		assert(myFile);
		// write out each block
		myFile.write(reinterpret_cast<char*>(m_filter), m_sizeInBytes);

		myFile.close();
		assert(myFile);
	}

	size_t getPop() const
	{
		size_t i, popBF = 0;
#pragma omp parallel for reduction(+ : popBF)
		for (i = 0; i < (m_size + 7) / 8; i++)
			popBF = popBF + popCnt(m_filter[i]);
		return popBF;
	}

	unsigned getHashNum() const { return m_hashNum; }

	unsigned getKmerSize() const { return m_kmerSize; }

	size_t getFilterSize() const { return m_size; }

	~BloomFilter() { delete[] m_filter; }
	// private:
	// BloomFilter(const BloomFilter& that); //to prevent copy construction

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
		m_filter = new unsigned char[m_sizeInBytes];
	}

	uint8_t* m_filter;
	size_t m_size;
	size_t m_sizeInBytes;
	unsigned m_hashNum;
	unsigned m_kmerSize;
};

#endif /* BLOOMFILTER_H_ */

//  Introducing member function pointer variables for multiHash:
vector<size_t> (BloomFilter::*fx1)(const char*) const = &BloomFilter::multiHash;
vector<size_t> (BloomFilter::*fx2)(const char*, uint64_t&, uint64_t&) const =
    &BloomFilter::multiHash;
vector<size_t> (BloomFilter::*fx3)(uint64_t&, uint64_t&, const char, const char) const =
    &BloomFilter::multiHash;

//  Introducing member function pointer variables for contains:
bool (BloomFilter::*fx4)(vector<size_t> const&) = &BloomFilter::contains;
bool (BloomFilter::*fx5)(const char*) = &BloomFilter::contains;
bool (BloomFilter::*fx6)(const char*, uint64_t&, uint64_t&) = &BloomFilter::contains;
bool (BloomFilter::*fx7)(uint64_t&, uint64_t&, const char, const char) = &BloomFilter::contains;

//  Introducing member function pointer variables for insert:
void (BloomFilter::*fx8)(vector<size_t> const&) = &BloomFilter::insert;
void (BloomFilter::*fx9)(const char*) = &BloomFilter::insert;
void (BloomFilter::*fx10)(const char*, uint64_t&, uint64_t&) = &BloomFilter::insert;
void (BloomFilter::*fx11)(uint64_t&, uint64_t&, const char, const char) = &BloomFilter::insert;

BOOST_PYTHON_MODULE(bloomfilter_ext)
{
	using namespace boost::python;

	class_<BloomFilter>("BloomFilter", init<size_t, unsigned, unsigned>())
	    .def(init<size_t, unsigned, unsigned, string>())
	    .def_readwrite("m_filter", &BloomFilter::m_filter)
	    .def_readwrite("m_size", &BloomFilter::m_size)
	    .def_readwrite("m_sizeInBytes", &BloomFilter::m_sizeInBytes)
	    .def_readwrite("m_hashNum", &BloomFilter::m_hashNum)
	    .def_readwrite("m_kmerSize", &BloomFilter::m_kmerSize)
	    .def("multiHash", fx1)
	    .def("multiHash", fx2)
	    .def("multiHash", fx3)
	    .def("contains", fx4)
	    .def("contains", fx5)
	    .def("contains", fx6)
	    .def("contains", fx7)
	    .def("insert", fx8)
	    .def("insert", fx9)
	    .def("insert", fx10)
	    .def("insert", fx11)
	    .def("storeFilter", &BloomFilter::storeFilter)
	    .def("getPop", &BloomFilter::getPop)
	    .def("getHashNum", &BloomFilter::getHashNum)
	    .def("getKmerSize", &BloomFilter::getKmerSize)
	    .def("getFilterSize", &BloomFilter::getFilterSize)
	    .def("initSize", &BloomFilter::initSize)

	    ;
}
