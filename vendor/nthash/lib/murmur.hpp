#ifndef MURMUR_HASH_H
#define MURMUR_HASH_H

#include <stdint.h>

//-----------------------------------------------------------------------------
// MurmurHash2, 64-bit versions, by Austin Appleby

// The same caveats as 32-bit MurmurHash2 apply here - beware of alignment
// and endian-ness issues if used across multiple platforms.

// 64-bit hash for 64-bit platforms

uint64_t MurmurHash64A ( const void * key, int len, unsigned int seed ) {
    const uint64_t m = 0xc6a4a7935bd1e995;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

    const uint64_t * data = (const uint64_t *)key;
    const uint64_t * end = data + (len/8);

    while (data != end)
    {
        uint64_t k = *data++;

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    const unsigned char * data2 = (const unsigned char*)data;

    switch(len & 7)
    {
    case 7:
        h ^= uint64_t(data2[6]) << 48; // fallthrough
    case 6:
        h ^= uint64_t(data2[5]) << 40; // fallthrough
    case 5:
        h ^= uint64_t(data2[4]) << 32; // fallthrough
    case 4:
        h ^= uint64_t(data2[3]) << 24; // fallthrough
    case 3:
        h ^= uint64_t(data2[2]) << 16; // fallthrough
    case 2:
        h ^= uint64_t(data2[1]) << 8; // fallthrough
    case 1:
        h ^= uint64_t(data2[0]);
        h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}


#endif
