# bloomfilter

The BTL C/C++ Common Bloom filters for bioinformatics projects, as well as any APIs created for other programming languages.

# usage example (C++)

Fast Bloom filter loading using the rolling hash function.

```C++
#include "BloomFilter.hpp"
#include <vector>
#include <string>
#include "vendor/ntHashIterator.hpp"

using namespace std;

int main(int argc, char** argv)
{
    /* test sequence */
    const string seq = "TAGAATCACCCAAAGA";
    /* k-mer size */
    const unsigned k = 5;
    /* number of Bloom filter hash functions */
    const unsigned numHashes = 4;
    /* size of Bloom filter (in bits) */
    const unsigned size = 1000;
	//Building the filter
	{
		/* init Bloom filter */
		BloomFilter bloom(size, numHashes, k);

		/* init rolling hash state and compute hash values for first k-mer */
		ntHashIterator itr(seq, numHashes, k);
		while (itr != itr.end()) {
			bloom.insert(*itr);
			++itr;
		}
		/* store the bloom filter */
		bloom.storeFilter("filterPathname.bf");
	}

	//After building
	{
		/* load the bloom filter */
		BloomFilter bloom("filterPathname.bf");

		/* query the bloom filter */

		/* init rolling hash state and compute hash values for first k-mer */
		ntHashIterator itr(seq, numHashes, k);
		while (itr != itr.end()) {
			bloom.contains(*itr);
			++itr;
		}
	}
	return 0;
}
```

Fast Counting Bloom filter loading using the rolling hash function.


```C++
#include "CountingBloomFilter.hpp"
#include <vector>
#include <string>
#include "vendor/ntHashIterator.hpp"

using namespace std;

int main(int argc, char** argv)
{
    /* test sequence */
    const string seq = "TAGAATCACCCAAAGA";
    /* k-mer size */
    const unsigned k = 5;
    /* number of Bloom filter hash functions */
    const unsigned numHashes = 4;
    /* size of Bloom filter (in bytes) */
    size_t size = 1000;
    /* counts to threshold bloom filter on*/
    const unsigned threshold = 1;
	//Building the filter
	{
		/* init Bloom filter */
		CountingBloomFilter<uint8_t> bloom(size, numHashes, k, threshold);

		/* init rolling hash state and compute hash values for first k-mer */
		ntHashIterator itr(seq, numHashes, k);
		while (itr != itr.end()) {
			bloom.insert(*itr);
			++itr;
		}
		/* store the bloom filter */
		bloom.storeFilter("filterPathname.bf");
	}

	//After building
	{
		/* load the bloom filter */
		CountingBloomFilter<uint8_t> bloom("filterPathname.bf", threshold);

		/* query the bloom filter */

		/* init rolling hash state and compute hash values for first k-mer */
		ntHashIterator itr(seq, numHashes, k);
		while (itr != itr.end()) {
			bloom.contains(*itr);
			++itr;
		}
	}
	return 0;
}
```

# files

* `BloomFilter.hpp`: main Bloom filter class
* `CountingBloomFilter.hpp`: Counting Bloom filter class
* `RollingHashIterator.h`: Enable rolling hashing on a string
* `RollingHash.h`: rolling hash interface (required by `RollingHashIterator.h`)
* `rolling.h`: rolling hash function (required by `BloomFilter.hpp` and `RollingHash.h`)
* `Tests/Unit`: unit tests
* `Tests/AdHoc`: ad-hoc tests

# unit tests

The unit tests may be compiled and run with:

	$ ./autogen.sh
	$ ./configure
	$ make check

To see more detailed output for the individual tests, run the binaries in `Tests/Unit` from the command line. (The ad-hoc tests in `Tests/AdHoc` may also be run in this way.)

# acknowledgements

This projects uses:
* [CATCH](https://github.com/philsquared/Catch) unit test framework for C/C++
* [nthash](https://github.com/bcgsc/ntHash) rolling hash implementation by Hamid Mohamadi
* [cpptoml](https://github.com/skystrife/cpptoml) TOML parser and serializer implemented by Chase Geigle

# Bloom filter file format

The specification of the Bloom filter file format is as follows:

1. magic header string
  * Description: bf magic string
  * Type: string
  * Value: BTLBloomFilter_v1 or BTLCountingBloomFilter_v1
2. header
  * Description: Plain header text
  * Type: string
  * Value:
    * size
      * Description: The size of Bloom filter
      * Type: size_t
      * Value:
    * sizeInBytes
      * Description: The size of Bloom filter in bytes
      * Type: size_t
      * Value:
    * hashNum
      * Description: number of hashes
      * Type: unsigned
      * Value:
    * kmerSize
      * Description: k-mer size
      * Type: unsigned
      * Value:
    * dFPR [optional]
      * Description: desired false positve rate
      * Type: double
      * Value:
    * nEntry [optional]
      * Description: number of entries
      * Type: uint64_t
      * Value:
    * tEntry [optional]
      * Description: total number of entries
      * Type: uint64_t
      * Value:
    * seed [optional] \(Not yet implimented\)
      * Description: initial seeds for different hashes
      * Type: uint64_t[nhash]
      * Value: [0,1, ..., nhash-1]
    * bitsPerCounter [optional]
      * Description: bits per each counter in the counting bloom filter
      * Type: unsigned
      * Value: 8
3. filter
  * Description: Bloom filter content
  * Type: uchar[sizeInBytes]
  * Value:
