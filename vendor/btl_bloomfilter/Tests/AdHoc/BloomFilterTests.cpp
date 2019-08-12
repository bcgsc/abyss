/*
 * BloomFilterTests.cpp
 * Unit Tests for hashmanager and bloomfilter classes
 *  Created on: Aug 14, 2012
 *      Author: cjustin
 */

#include "BloomFilter.hpp"
#include "vendor/ntHashIterator.hpp"
#include <assert.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#if _OPENMP
#include <omp.h>
#endif

using namespace std;

// returns memory of program in kb
unsigned
memory_usage()
{
	unsigned mem = 0;
	ifstream proc("/proc/self/status");
	string s;
	while (getline(proc, s), !proc.fail()) {
		if (s.substr(0, 6) == "VmSize") {
			stringstream convert(s.substr(s.find_last_of('\t'), s.find_last_of('k') - 1));
			if (!(convert >> mem)) {
				return 0;
			}
			return mem;
		}
	}
	return mem;
}

int
main()
{
	// memory usage from before
	int memUsage = memory_usage();

	size_t filterSize = 1000000;

	const unsigned numHashes = 3;
	const unsigned k = 4;
	const char* seq = "ACGTAC";

	BloomFilter filter(filterSize, numHashes, k);

	// insert k-mers ACGT, CGTA, GTAC

	ntHashIterator insertIt(seq, numHashes, k);
	while (insertIt != insertIt.end()) {
		filter.insert(*insertIt);
		++insertIt;
	}

	// check that k-mers were correctly inserted

	ntHashIterator queryIt(seq, numHashes, k);
	while (queryIt != queryIt.end()) {
		assert(filter.contains(*queryIt));
		++queryIt;
	}

	// should be size of bf (amortized)
	cout << memory_usage() - memUsage << "kb" << endl;

	cout << "de novo bf tests done" << endl;

	// Check storage can occur properly
	string filename = "/tmp/bloomFilter.bf";
	filter.storeFilter(filename);
	ifstream ifile(filename.c_str());
	assert(ifile.is_open());
	// Read header line by line
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
	// Get header and file size
	size_t currPos = ifile.tellg();
	ifile.seekg(0, ios::end);        // move to end of file
	size_t fileSize = ifile.tellg(); // file size in bytes
	// file size - header size should be same as filter size (Round to block size)
	if (filterSize % 64 > 0) {
		assert((filterSize + (64 - (filterSize % 64))) == (fileSize - currPos) * 8);
	} else {
		assert(filterSize == (fileSize - currPos) * 8);
	}
	ifile.close();

	// should be roughly same size still (amortized)
	cout << memory_usage() - memUsage << "kb" << endl;

	// check loading of stored filter
	BloomFilter filter2(filename);

	// should be double size of bf (amortized)
	cout << memory_usage() - memUsage << "kb" << endl;

	// Check if loaded filter is able to report expected results

	ntHashIterator queryIt2(seq, numHashes, k);
	while (queryIt2 != queryIt2.end()) {
		assert(filter2.contains(*queryIt2));
		++queryIt2;
	}

	cout << "premade bf tests done" << endl;

#if __linux__

	// memory leak tests
	BloomFilter* filter3 = new BloomFilter(filterSize, 5, 20);

	size_t tempMem = memory_usage() - memUsage;

	cout << memory_usage() - memUsage << "kb" << endl;
	delete (filter3);
	cout << memory_usage() - memUsage << "kb" << endl;
	assert(tempMem != memory_usage() - memUsage);

	//	vector<vector<int> > test(1, vector<int>(1, 1));

	//	vector<BloomFilter> test;
	//	test.push_back(filter2);

	cout << "memory leak prevention tests done" << endl;
	cout << memory_usage() - memUsage << "kb" << endl;

#endif

	remove(filename.c_str());

	//	//check parallelized code speed
	//	cout << "testing code parallelization" << endl;
	//	double start_s = omp_get_wtime();
	//
	//	cout << "parallelization" << endl;
	//	for( int i =0; i < 10000000; i++)
	//	{
	//		vector<size_t> values = multiHash("ATCGGGTCATCAACCAATAA", 5, 20);
	//		filter.contains(values);
	//	}
	//	double stop_s=omp_get_wtime();
	//	cout << "time: " << stop_s-start_s << endl;
	//
	//	start_s = omp_get_wtime();
	//	cout << "non parallelization" << endl;
	//	for( int i =0; i < 10000000; i++)
	//	{
	//		vector<size_t> values = multiHashNonPara("ATCGGGTCATCAACCAATAA", 5, 20);
	//		filter.contains(values);
	//	}
	//	stop_s=omp_get_wtime();
	//	cout << "time: " << stop_s-start_s << endl;

	cout << "done" << endl;
	return 0;
}
