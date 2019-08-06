%module BloomFilter
%include "std_string.i"
%include "stdint.i"
%include "std_vector.i"
/*%apply const uint64_t& { uint64_t & };*/
/*%apply uint64_t& INOUT {uint64_t& fhVal, uint64_t& rhVal};*/
namespace std {
   %template(SizetVector) vector<size_t>;
}

%{
#include "../KmerBloomFilter.hpp"
#include "../vendor/ntHashIterator.hpp"
#include "../BloomFilterUtil.h"
%}

%rename(BloomFilter) KmerBloomFilter;

using namespace std;

class KmerBloomFilter {
public:
        KmerBloomFilter();
        ~KmerBloomFilter();
        KmerBloomFilter(uint64_t filterSize, unsigned hashNum, unsigned kmerSize);
        KmerBloomFilter(const string &filterFilePath);

        void insert(vector<uint64_t> const &precomputed);
        void insert(const char* kmer);

        bool contains(vector<uint64_t> const &values);
        bool contains(const char* kmer);

        void storeFilter(string const &filterFilePath);
        uint64_t getPop();
        unsigned getHashNum();
        unsigned getKmerSize();
        uint64_t getFilterSize();
};

/*
class ntHashIterator {
public:
    ntHashIterator();
    ~ntHashIterator();
    ntHashIterator(const string& seq, unsigned numHashes, unsigned k);
 
    const uint64_t* operator*();
    const uint64_t* operator->();

    bool operator==(const ntHashIterator& it);
    bool operator!=(const ntHashIterator& it);
    
    ntHashIterator& operator++();
    static const ntHashIterator end();
};
*/

void insertSeq(KmerBloomFilter &bloom, const string& seq, unsigned numHashes, unsigned k);
