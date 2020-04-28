
/**
 * Konnector Bloom Filter compatability layer
 * Copyright 2020 Johnahtan Wong
 */
#ifndef KONNECTORBLOOMFILTER_H
#define KONNECTORBLOOMFILTER_H 1

#include <execinfo.h>

#include "Bloom/Bloom.h"
#include "BloomDBG/RollingHashIterator.h"
#include "Common/BitUtil.h"
#include "Common/IOUtil.h"
#include "Common/Kmer.h"
#include "vendor/btl_bloomfilter/BloomFilter.hpp"
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <vector>

#define BT_BUF_SIZE 100

/** A Bloom filter. */
class KonnectorBloomFilter : public BloomFilter
{
  public:
	/** Constructor. */
	KonnectorBloomFilter()
	  : BloomFilter()
	{}

	/** Constructor. */
	KonnectorBloomFilter(size_t n, unsigned k)
	  : BloomFilter{ n, 1, k }
	{}

	~KonnectorBloomFilter() {}

	/** Return the size of the bit array. */
	size_t size() const { return getFilterSize(); }

	/** Return the population count, i.e. the number of set bits. */
	size_t popcount() const { return getPop(); }

	/** Return the estimated false positive rate */
	double FPR() const { return getFPR(); }

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	bool contains(vector<size_t> const& precomputed) const
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
	bool contains(const size_t precomputed[]) const
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

	/** Return whether the specified bit is set. */
	bool operator[](const size_t precomputed[]) const { return contains(precomputed); }

	/** Return whether the specified bit is set. */
	bool operator[](size_t i) const
	{
		size_t foo[1] = { i };
		return contains(foo);
	}

	/** Return whether the object is present in this set. */
	bool operator[](const Bloom::key_type& key) const
	{
		RollingHashIterator it(key.str().c_str(), 1, key.length());
		return contains(*it);
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	void insert(const size_t precomputed[])
	{

		// iterates through hashed values adding it to the filter
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t normalizedValue = precomputed[i] % m_size;
			__sync_or_and_fetch(
			    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
		}
	}

	void insert(vector<size_t> const& precomputed)
	{

		// iterates through hashed values adding it to the filter
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t normalizedValue = precomputed.at(i) % m_size;
			__sync_or_and_fetch(
			    &m_filter[normalizedValue / bitsPerChar], bitMask[normalizedValue % bitsPerChar]);
		}
	}

	/** Add the object with the specified index to this set. */
	void insert(size_t i)
	{
		size_t foo[1] = { i };
		insert(foo);
	}

	/** Add the object to this set. */
	void insert(const Bloom::key_type& key)
	{
		RollingHashIterator it(key.str().c_str(), 1, key.length());
		insert(*it);
	}

	/** Operator for reading a bloom filter from a stream. */
	friend std::istream& operator>>(std::istream& in, KonnectorBloomFilter& o)
	{
		o.read(in, BITWISE_OVERWRITE);
		return in;
	}

	/** Operator for writing the bloom filter to a stream. */
	friend std::ostream& operator<<(std::ostream& out, const KonnectorBloomFilter& o)
	{
		o.writeHeader(out);
		// NOLINTNEXTLINE(google-readability-casting)
		out.write(reinterpret_cast<char*>(o.m_filter), o.m_sizeInBytes);
		return out;
	}

	void loadHeaderKonnector(std::istream& file)
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
		m_incomingSize = *bloomFilterTable->get_as<size_t>("BloomFilterSize");
		m_hashNum = *bloomFilterTable->get_as<unsigned>("HashNum");
		m_kmerSize = *bloomFilterTable->get_as<unsigned>("KmerSize");
		m_sizeInBytes = *bloomFilterTable->get_as<size_t>("BloomFilterSizeInBytes");
		m_dFPR = *bloomFilterTable->get_as<double>("dFPR");
		m_nEntry = *bloomFilterTable->get_as<uint64_t>("nEntry");
		m_tEntry = *bloomFilterTable->get_as<uint64_t>("Entry");
	}

	/** Read a bloom filter from a stream. */
	void read(std::istream& in, BitwiseOp readOp = BITWISE_OVERWRITE)
	{
		loadHeaderKonnector(in);

		if (m_size != m_incomingSize) {
			if (readOp == BITWISE_OVERWRITE) {
				resize(m_incomingSize);
			} else {
				std::cerr << "error: can't union/intersect bloom filters with "
				          << "different sizes\n";
				exit(EXIT_FAILURE);
			}
		}

		readBits(in, reinterpret_cast<char*>(m_filter), m_size, 0, readOp);
		assert(in);
	}

	/** Write a bloom filter to a stream. */
	void write(std::ostream& out) const
	{
		writeHeader(out);
		out.write(reinterpret_cast<char*>(m_filter), m_sizeInBytes);
		assert(out);
	}

	/** Resize the bloom filter (wipes the current data) */
	void resize(size_t size)
	{
		initSize(size);
		m_size = size;
	}

  protected:
	size_t m_incomingSize;
};

#endif
