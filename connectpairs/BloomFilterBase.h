#ifndef BLOOMFILTERBASE_H
#define BLOOMFILTERBASE_H

#include "Common/Kmer.h"
#include "Common/HashFunction.h"
#include "Common/Uncompress.h"
#include "DataLayer/FastaReader.h"
#include <iostream>

class BloomFilterBase
{
public:

	static const unsigned BLOOM_VERSION = 1;
	typedef Kmer key_type;

	/** Constructor. */
	BloomFilterBase() { }

	/** Add the object to this set. */
	virtual void insert(const key_type& key) = 0;

	/** Return whether the object is present in this set. */
	virtual bool operator[](const key_type& key) const = 0;

	/** Return the size of the bit array. */
	virtual size_t size() const = 0;

	/** Return the population count, i.e. the number of set bits. */
	virtual size_t popcount() const = 0;

	/** Return the estimated false positive rate */
	virtual double FPR() const = 0;

	/** Return the hash value of this object. */
	static size_t hash(const key_type& key)
	{
		return hashmem(&key, sizeof key);
	}

	void loadSeq(unsigned k, const std::string& seq)
	{
		for (size_t i = 0; i < seq.size() - k + 1; ++i) {
			std::string kmer = seq.substr(i, k);
			size_t pos = kmer.find_last_not_of("ACGTacgt");
			if (pos == std::string::npos) {
				Kmer u(kmer);
				insert(u);
				u.reverseComplement();
				insert(u);
			} else
				i += pos;
		}
	}

	void loadFile(unsigned k, const std::string& path, bool verbose = false)
	{
		assert(!path.empty());
		if (verbose)
			std::cerr << "Reading `" << path << "'...\n";
		FastaReader in(path.c_str(), FastaReader::NO_FOLD_CASE);
		uint64_t count = 0;
		for (std::string seq; in >> seq; count++) {
			if (verbose && count % LOAD_PROGRESS_STEP == 0)
				std::cerr << "Loaded " << count << " reads into bloom filter\n";
			loadSeq(k, seq);
		}
		assert(in.eof());
		if (verbose)
			std::cerr << "Loaded " << count << " reads into bloom filter\n";
	}

private:

	static const unsigned LOAD_PROGRESS_STEP = 100000;

};

#endif
