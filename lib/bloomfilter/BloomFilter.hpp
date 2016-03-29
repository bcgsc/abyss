/*
 *
 * BloomFilter.hpp
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_
#include <string>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include "rolling.h"

using namespace std;

static const uint8_t bitsPerChar = 0x08;
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
    0x40, 0x80 };

inline unsigned popCnt(unsigned char x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
    & 0xf;
}

class BloomFilter {
public:

    /*
     * Default constructor.
    */
    BloomFilter() : m_filter(NULL), m_size(0), m_sizeInBytes(0),
        m_hashNum(0), m_kmerSize(0) {}

    /* De novo filter constructor.
     *
     * preconditions:
     * filterSize must be a multiple of 64
     * kmerSize refers to the number of bases the kmer has
     * k-mers supplied to this object should be binary (2 bits per base)
     */
    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize) :
    m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        initSize(m_size);
        memset(m_filter, 0, m_sizeInBytes);
    }

    /*
     * Loads the filter (file is a .bf file) from path specified
     */
    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize,
                string const &filterFilePath) :
    m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        initSize(m_size);

        FILE *file = fopen(filterFilePath.c_str(), "rb");
        if (file == NULL) {
            cerr << "file \"" << filterFilePath << "\" could not be read."
            << endl;
            exit(1);
        }

        long int lCurPos = ftell(file);
        fseek(file, 0, 2);
        size_t fileSize = ftell(file);
        fseek(file, lCurPos, 0);
        if (fileSize != m_sizeInBytes) {
            cerr << "Error: " << filterFilePath
            << " does not match size given by its information file. Size: "
            << fileSize << " vs " << m_sizeInBytes << " bytes." << endl;
            exit(1);
        }

        size_t countRead = fread(m_filter, fileSize, 1, file);
        if (countRead != 1 && fclose(file) != 0) {
            cerr << "file \"" << filterFilePath << "\" could not be read."
            << endl;
            exit(1);
        }
    }

    /*
     * Accepts a list of precomputed hash values. Faster than rehashing each time.
     */
    void insert(vector<size_t> const &precomputed) {

        //iterates through hashed values adding it to the filter
        for (size_t i = 0; i < m_hashNum; ++i) {
            size_t normalizedValue = precomputed.at(i) % m_size;
            __sync_or_and_fetch(&m_filter[normalizedValue / bitsPerChar],
                                bitMask[normalizedValue % bitsPerChar]);
        }
    }

    void insert(const char* kmer) {
        uint64_t hVal = getChval(kmer, m_kmerSize);
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t normalizedValue = (rol(varSeed, i) ^ hVal) % m_size;
            __sync_or_and_fetch(&m_filter[normalizedValue / bitsPerChar],
                                bitMask[normalizedValue % bitsPerChar]);
        }
    }

    /*
     * Accepts a list of precomputed hash values. Faster than rehashing each time.
     */
    bool contains(vector<size_t> const &values) const {
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
    bool contains(const char* kmer) const {
        uint64_t hVal = getChval(kmer, m_kmerSize);
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
    void storeFilter(string const &filterFilePath) const {
        ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

        cerr << "Storing filter. Filter is " << m_sizeInBytes << "bytes."
        << endl;

        assert(myFile);
        //write out each block
        myFile.write(reinterpret_cast<char*>(m_filter), m_sizeInBytes);

        myFile.close();
        assert(myFile);
    }

    size_t getPop() const {
        size_t i, popBF=0;
#pragma omp parallel for reduction(+:popBF)
        for(i=0; i<(m_size + 7)/8; i++)
            popBF = popBF + popCnt(m_filter[i]);
        return popBF;
    }

    unsigned getHashNum() const {
        return m_hashNum;
    }

    unsigned getKmerSize() const {
        return m_kmerSize;
    }

    size_t getFilterSize() const { return m_size; }

    ~BloomFilter() {
        delete[] m_filter;
    }
private:
    BloomFilter(const BloomFilter& that); //to prevent copy construction

    /*
     * Checks filter size and initializes filter
     */
    void initSize(size_t size) {
        if (size % 8 != 0) {
            cerr << "ERROR: Filter Size \"" << size
            << "\" is not a multiple of 8." << endl;
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
