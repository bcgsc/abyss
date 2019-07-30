#include <iomanip>
#include <iostream>
#include <string>

#include "BloomFilter.hpp"
#include "vendor/ntHashIterator.hpp"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdint.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace opt {
/** default size of the Bloom filter in bits (1MB) */
size_t bloomBits = 1024 * 1024 * 8;
unsigned kmerLen = 64;
unsigned ibits = 64;
unsigned nhash = 5;
}

using namespace std;

static const unsigned char b2r[256] = { 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 0
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 1
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 2
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 3
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', // 4   'A' 'C' 'G'
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', // 5   'T'
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', // 6   'a' 'c' 'g'
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', // 7   't'
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 8
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 9
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 10
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 11
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 12
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 13
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 14
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 15
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N' };

void
getCanon(std::string& bMer)
{
	int p = 0, hLen = (opt::kmerLen - 1) / 2;
	while (bMer[p] == b2r[(unsigned char)bMer[opt::kmerLen - 1 - p]]) {
		++p;
		if (p >= hLen)
			break;
	}
	if (bMer[p] > b2r[(unsigned char)bMer[opt::kmerLen - 1 - p]]) {
		for (int lIndex = p, rIndex = opt::kmerLen - 1 - p; lIndex <= rIndex; ++lIndex, --rIndex) {
			char tmp = b2r[(unsigned char)bMer[rIndex]];
			bMer[rIndex] = b2r[(unsigned char)bMer[lIndex]];
			bMer[lIndex] = tmp;
		}
	}
}

void
loadSeq(BloomFilter& BloomFilterFilter, const string& seq)
{
	if (seq.size() < opt::kmerLen)
		return;

	ntHashIterator insertIt(seq, BloomFilterFilter.getHashNum(), BloomFilterFilter.getKmerSize());
	while (insertIt != insertIt.end()) {
		BloomFilterFilter.insert(*insertIt);
		++insertIt;
	}
}

void
loadSeqr(BloomFilter& BloomFilterFilter, const string& seq)
{
	if (seq.size() < opt::kmerLen)
		return;
	string kmer = seq.substr(0, opt::kmerLen);
	ntHashIterator itr(seq, opt::kmerLen, opt::nhash);
	while (itr != itr.end()) {
		BloomFilterFilter.insert(*itr);
		++itr;
	}
}

void
loadBf(BloomFilter& BloomFilterFilter, const char* faqFile)
{
	ifstream uFile(faqFile);
	bool good = true;
#pragma omp parallel
	for (string line, hline; good;) {
#pragma omp critical(uFile)
		{
			good = static_cast<bool>(getline(uFile, hline));
			good = static_cast<bool>(getline(uFile, line));
			// good = getline(uFile, hline);
			// good = getline(uFile, hline);
		}
		if (good)
			loadSeqr(BloomFilterFilter, line);
	}
	uFile.close();
}

int
main(int argc, const char* argv[])
{
	/*BloomFilter myFilter(40857600000, 2, 30);
	string mystr="AGAGACGTGCATCGGGTCATCAACCAATAT";
	myFilter.insert(mystr.c_str());
	if (myFilter.search(mystr.c_str()))
	    cerr << mystr << " is in Bloom filter.\n";
	else
	    cerr << mystr << " is not in Bloom filter.\n";
	return 0;*/
	if (argc < 2)
		cerr << "error!\n";
#ifdef _OPENMP
	double sTime = omp_get_wtime();
#endif
	/*
	 * Note: The previous Bloom filter size here was
	 * 48,857,600,000 bits (~6GB).  I reduced it to 1024*1024*8
	 * (1MB) by default, so that it would run safely on any machine.
	 * -BV
	 */
	BloomFilter myFilter(opt::bloomBits, opt::nhash, opt::kmerLen);
	loadBf(myFilter, argv[1]);
	cerr << "|popBF|=" << myFilter.getPop() << " ";
#ifdef _OPENMP
	cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
#else
	cerr << "\n";
#endif
	// myFilter.store("filter3.bf");

	/*BloomFilter filter2(40857600000, 5, 30, "filter1.bf");
	cerr << "|popBF|=" << filter2.getPop() << " ";
	cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";

	filter2.store("filter2.bf");*/

	return 0;
}
