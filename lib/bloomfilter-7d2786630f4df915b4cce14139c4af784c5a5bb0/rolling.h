#ifndef ROLLING_HASH_H
#define ROLLING_HASH_H

#include <stdint.h>

// offset for the complement base in the random seeds table
const int cpOff = -20;

// seed for gerenerating multiple hash values
const uint64_t varSeed = 2577914034309095328ul;

// 64-bit random seed table corresponding to bases and their complements
static const uint64_t seedTab[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, // 0..7
    0, 0, 0, 0, 0, 0, 0, 0, // 8..15
    0, 0, 0, 0, 0, 0, 0, 0, // 16..23
    0, 0, 0, 0, 0, 0, 0, 0, // 24..31
    0, 0, 0, 0, 0, 0, 0, 0, // 32..39
    0, 0, 0, 0, 0, 2978368046464386134ul, 0, 2319985823310095140ul, // 40..47
    0, 0, 0, 3572411708064410444ul, 0, 0, 0, 0, // 48..55
    0, 0, 0, 0, 0, 0, 0, 0, // 56..63
    4362857412768957556ul, 4362857412768957556ul, 0, 3572411708064410444ul, 0, 0, 0, 2319985823310095140ul, // 64..71
    0, 0, 0, 0, 0, 2978368046464386134ul, 0, 2319985823310095140ul, // 72..79
    0, 0, 0, 3572411708064410444ul, 2978368046464386134ul, 0, 0, 0, // 80..87
    0, 0, 0, 0, 0, 0, 0, 0, // 88..95
    4362857412768957556ul, 4362857412768957556ul, 0, 3572411708064410444ul, 0, 0, 0, 2319985823310095140ul, // 96..103
    0, 0, 0, 0, 0, 0, 0, 0, // 104..111
    0, 0, 0, 0, 2978368046464386134ul, 0, 0, 0, // 112..119
    0, 0, 0, 0, 0, 0, 0, 0, // 120..127
    0, 0, 0, 0, 0, 0, 0, 0, // 128..135
    0, 0, 0, 0, 0, 0, 0, 0, // 136..143
    0, 0, 0, 0, 0, 0, 0, 0, // 144..151
    0, 0, 0, 0, 0, 0, 0, 0, // 152..159
    0, 0, 0, 0, 0, 0, 0, 0, // 160..167
    0, 0, 0, 0, 0, 0, 0, 0, // 168..175
    0, 0, 0, 0, 0, 0, 0, 0, // 176..183
    0, 0, 0, 0, 0, 0, 0, 0, // 184..191
    0, 0, 0, 0, 0, 0, 0, 0, // 192..199
    0, 0, 0, 0, 0, 0, 0, 0, // 200..207
    0, 0, 0, 0, 0, 0, 0, 0, // 208..215
    0, 0, 0, 0, 0, 0, 0, 0, // 216..223
    0, 0, 0, 0, 0, 0, 0, 0, // 224..231
    0, 0, 0, 0, 0, 0, 0, 0, // 232..239
    0, 0, 0, 0, 0, 0, 0, 0, // 240..247
    0, 0, 0, 0, 0, 0, 0, 0  // 248..255
};

// rotate "v" to the left by "s" positions
inline uint64_t rol(const uint64_t v, const int s) {
    return (v << s) | (v >> (64 - s));
}

// rotate "v" to the right by "s" positions
inline uint64_t ror(const uint64_t v, const int s) {
    return (v >> s) | (v << (64 - s));
}

// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
inline uint64_t getFhval(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    return hVal;
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t getRhval(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]+cpOff], i);
    return hVal;
}

// cannonical hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t getChval(const char * kmerSeq, const unsigned k) {
    uint64_t fhVal = getFhval(kmerSeq, k);
    uint64_t rhVal = getRhval(kmerSeq, k);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// initialize forward-strand hash value of the first kmer, i.e. fhval(kmer_0)
inline uint64_t initHashes(const char * kmerSeq, const unsigned k) {
    return getFhval(kmerSeq, k);
}

// initialize cannonical hash value of the first kmer, i.e. chval(kmer_0)
inline uint64_t initHashes(const char * kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = getFhval(kmerSeq, k);
    rhVal = getRhval(kmerSeq, k);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// recursive forward-strand hash value for next k-mer
inline uint64_t rollHashesRight(const uint64_t fhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    return(rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn]);
}

// recursive cannonical hash value for next k-mer
inline uint64_t rollHashesRight(uint64_t& fhVal, uint64_t& rhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    fhVal = rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn];
    rhVal = ror(rhVal, 1) ^ ror(seedTab[charOut+cpOff], 1) ^ rol(seedTab[charIn+cpOff], k-1);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// recursive forward-strand hash value for prev k-mer
inline uint64_t rollHashesLeft(const uint64_t fhVal, const unsigned char charIn, const unsigned char charOut, const unsigned k) {
    return(ror(fhVal, 1) ^ ror(seedTab[charOut], 1) ^ rol(seedTab[charIn], k-1));
}

// recursive canonical hash value for prev k-mer
inline uint64_t rollHashesLeft(uint64_t& fhVal, uint64_t& rhVal, const unsigned char charIn, const unsigned char charOut, const unsigned k) {
    fhVal = ror(fhVal, 1) ^ ror(seedTab[charOut], 1) ^ rol(seedTab[charIn], k-1);
    rhVal = rol(rhVal, 1) ^ rol(seedTab[charOut+cpOff], k) ^ seedTab[charIn+cpOff];
    return (rhVal<fhVal)? rhVal : fhVal;
}

#endif
