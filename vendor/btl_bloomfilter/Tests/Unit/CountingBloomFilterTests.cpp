/*
 * CountingBloomFilterTests.cpp
 * Unit Tests for CountingBloomFilter class
 * Adapted from BloomFilterTests.cpp
 *  Created on: July 15, 2019
 *      Author: Johnathan Wong
 */

/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN

/* lightweight unit test framework */
#include "CountingBloomFilter.hpp"
#include "vendor/catch.hpp"
#include "vendor/ntHashIterator.hpp"

#include <assert.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>

using namespace std;

/** create a uniquely-named temp file (tedious!) */
string
createTempFile()
{
	const unsigned MAX_FILENAME_SIZE = 1024;
	const char* fileTemplate = "/XXXXXX";
	char filename[MAX_FILENAME_SIZE + 1];

	/* allow override of default tmp dir */
	char* tmpdir = getenv("TMPDIR");
	if (tmpdir)
		strcpy(filename, tmpdir);
	else
		strcpy(filename, "/tmp");

	assert(strlen(filename) + strlen(fileTemplate) <= MAX_FILENAME_SIZE);
	strcat(filename, fileTemplate);

	int fd = mkstemp(filename);
	if (fd == -1) {
		perror("failed to create temp file");
		exit(EXIT_FAILURE);
	}
	close(fd);

	return string(filename);
}

TEST_CASE("test fixture", "[CountingBloomFilter]")
{
	/*
	 * NOTES:
	 * - The SECTION blocks below are separate tests that share the
	 * same setup code
	 * - The common setup code is _re-run_ before each SECTION
	 * - In unit test terminology, this type of setup is known as a
	 * "test fixture"
	 * - See
	 * https://github.com/philsquared/Catch/blob/master/docs/tutorial.md#test-cases-and-sections for
	 * details
	 */

	/* START COMMON SETUP CODE */

	const size_t filterSize = 100001;
	const unsigned numHashes = 5;
	const unsigned threshold = 1;
	const unsigned k = 8;
	const char* seq = "ACGTACACTGGACTGAGTCT";

	CountingBloomFilter<uint8_t> filter(filterSize, numHashes, k, threshold);

	/* insert k-mers ACGT, CGTA, GTAC */

	ntHashIterator insertIt(seq, numHashes, k);
	while (insertIt != insertIt.end()) {
		filter.insert(*insertIt);
		++insertIt;
	}

	CountingBloomFilter<uint64_t> filter_64bit(filterSize, numHashes, k, threshold);

	/* insert k-mers ACGT, CGTA, GTAC */

	ntHashIterator insertIt2(seq, numHashes, k);
	while (insertIt2 != insertIt2.end()) {
		filter_64bit.insert(*insertIt2);
		++insertIt2;
	}

	/* END COMMON SETUP CODE */

	SECTION("query elements")
	{
		/* check that k-mers were correctly inserted */

		ntHashIterator queryIt(seq, numHashes, k);
		while (queryIt != queryIt.end()) {
			assert(filter.contains(*queryIt));
			assert(filter_64bit.contains(*queryIt));
			++queryIt;
		}

		/* check that k-mers were not incorrectly inserted */

		string seq2;
		string DNA = "ATCG";
		srand(time(0));
		for (int i = 0; i < 60; i++) {
			seq2 += DNA[rand() % 4];
		}
		ntHashIterator queryIt2(seq2.c_str(), numHashes, k);
		while (queryIt2 != queryIt2.end()) {
			assert(!filter.contains(*queryIt2));
			assert(!filter_64bit.contains(*queryIt2));
			++queryIt2;
		}
	}

	SECTION("save/load 8 bit Bloom file")
	{
		/* write filter */

		string filename = createTempFile();
		filter.storeFilter(filename);
		ifstream ifile(filename.c_str());

		/* check size of newly-created file */

		assert(ifile.is_open());

		/* File header has no fixed element but "[HeaderEnd]"
		   so check for "[HeaderEnd]" in file*/

		std::string headerEnd = "[HeaderEnd]";
		std::string line;
		bool headerEndCheck = false;
		while (std::getline(ifile, line)) {
			if (line == headerEnd) {
				headerEndCheck = true;
				break;
			}
		}
		assert(headerEndCheck);

		size_t currPos = ifile.tellg();
		ifile.seekg(0, ios::end);
		size_t endPos = ifile.tellg();

		// file size - header should be same as filter size rounded to multiple of 64

		int remainder = filterSize % 8;
		assert((endPos - currPos) == (filterSize + 8 - remainder));
		ifile.close();

		/* check loading of stored filter */

		CountingBloomFilter<uint8_t> filter2(filename, threshold);

		// Checking if sizeInBytes correspond to filesize - header

		assert((endPos - currPos) == filter2.sizeInBytes());

		// In a 8bit filter, size should be equal to size_in_bytes

		assert(filter2.size() == filter2.sizeInBytes());

		/* check if loaded filter is able to report expected results */

		ntHashIterator queryIt(seq, numHashes, k);
		while (queryIt != queryIt.end()) {
			assert(filter2.contains(*queryIt));
			++queryIt;
		}

		/* cleanup */

		remove(filename.c_str());
	}

	SECTION("save/load 64 bit Bloom file")
	{
		/* write filter */

		string filename = createTempFile();
		filter_64bit.storeFilter(filename);
		ifstream ifile(filename.c_str());

		/* check size of newly-created file */

		assert(ifile.is_open());
		/* File header has no fixed element but "[HeaderEnd]"
		   so check for "[HeaderEnd]" in file*/
		std::string headerEnd = "[HeaderEnd]";
		std::string line;
		bool headerEndCheck = false;
		while (std::getline(ifile, line)) {
			if (line == headerEnd) {
				headerEndCheck = true;
				break;
			}
		}
		assert(headerEndCheck);

		size_t currPos = ifile.tellg();
		ifile.seekg(0, ios::end);
		size_t endPos = ifile.tellg();

		// file size - header should be same as filter size rounded to multiple of 64

		int remainder = filterSize % 8;
		assert((endPos - currPos) == (filterSize + 8 - remainder));

		ifile.close();

		/* check loading of stored filter */

		CountingBloomFilter<uint64_t> filter_64bit2(filename, threshold);

		// Checking if sizeInBytes correspond to filesize - header

		assert((endPos - currPos) == filter_64bit2.sizeInBytes());

		// In a 64bit filter, size * 8 should be equal to size_in_bytes

		assert((filter_64bit2.size() * 8) == filter_64bit2.sizeInBytes());

		/* check if loaded filter is able to report expected results */

		ntHashIterator queryIt(seq, numHashes, k);
		while (queryIt != queryIt.end()) {
			assert(filter_64bit2.contains(*queryIt));
			++queryIt;
		}

		/* cleanup */

		remove(filename.c_str());
	}

} /* end test fixture */
