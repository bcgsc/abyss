#include "BloomFilter.hpp"
#include <string>
#include <vector>

using namespace std;

/** stores state between calls to rolling hash */
struct RollingHashState
{
	/* seed hash value for current k-mer */
	uint64_t hash;
	/* seed hash value for reverse complement of current k-mer */
	uint64_t rcHash;
};

int
main(int argc, char** argv)
{
	/* test sequence */
	const string seq = "TAGAATCACCCAAAGA";
	/* k-mer size */
	const unsigned k = 5;
	/* number of Bloom filter hash functions */
	const unsigned numHashes = 4;
	/* size of Bloom filter (in bits) */
	const unsigned size = 1000;
	/* hash values for current k-mer */
	vector<size_t> hashes;

	/* init Bloom filter */
	BloomFilter bloom(size, numHashes, k);

	/* init rolling hash state and compute hash values for first k-mer */
	RollingHashState state;
	string kmer0 = seq.substr(0, k);
	hashes = bloom.multiHash(kmer0.c_str(), state.hash, state.rcHash);

	/* load k-mers into Bloom filter using rolling hash */
	for (unsigned i = 1; i < seq.length() - k + 1; ++i) {
		/* "roll" hash values right to current k-mer */
		char charOut = seq[i - 1];
		char charIn = seq[i + k - 1];
		hashes = bloom.multiHash(state.hash, state.rcHash, charOut, charIn);
		/* insert current k-mer into Bloom filter */
		bloom.insert(hashes);
	}

	return 0;
}
