/**
 * A cascading Bloom filter
 * Copyright 2013 Shaun Jackman
 */
#ifndef CascadingBLOOMFILTER_H
#define CascadingBLOOMFILTER_H 1

#include "Bloom/Bloom.h"
#include "BloomFilter.h"
#include <vector>

/** A Cascading Bloom filter. */
class CascadingBloomFilter : BloomFilter
{
  public:

	/** The maximum count of an element in this multiset. */
	static const unsigned MAX_COUNT = 2;

	/** Constructor */
	CascadingBloomFilter() {}

	/** Constructor */
	CascadingBloomFilter(size_t n)
	{
		for (unsigned i = 0; i < MAX_COUNT; i++)
			m_data.push_back(new BloomFilter(n));
	}

	/** Destructor */
	~CascadingBloomFilter()
	{
		typedef std::vector<BloomFilter*>::iterator Iterator;
		for (Iterator i = m_data.begin(); i != m_data.end(); i++) {
			assert(*i != NULL);
			delete *i;
		}
	}

	/** Return the size of the bit array. */
	size_t size() const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->size();
	}

	/** Return the number of elements with count >= MAX_COUNT. */
	size_t popcount() const
	{
		assert(m_data.back() != NULL);
		return m_data.back()->popcount();
	}

	/** Return the estimated false positive rate */
	double FPR() const
	{
		return (double)popcount() / size();
	}

	/** Return whether the element with this index has count >=
	 * MAX_COUNT.
	 */
	bool operator[](size_t i) const
	{
		assert(m_data.back() != NULL);
		return (*m_data.back())[i];
	}

	/** Return whether this element has count >= MAX_COUNT. */
	bool operator[](const Bloom::key_type& key) const
	{
		assert(m_data.back() != NULL);
		return (*m_data.back())[Bloom::hash(key) % m_data.back()->size()];
	}

	/** Add the object with the specified index to this multiset. */
	void insert(size_t index)
	{
		for (unsigned i = 0; i < MAX_COUNT; ++i) {
			assert(m_data.at(i) != NULL);
			if (!(*m_data[i])[index]) {
				m_data[i]->insert(index);
				break;
			}
		}
	}

	/** Add the object to this Cascading multiset. */
	void insert(const Bloom::key_type& key)
	{
		assert(m_data.back() != NULL);
		insert(Bloom::hash(key) % m_data.back()->size());
	}

	/** Get the Bloom filter for a given level */
	BloomFilter& getBloomFilter(unsigned level)
	{
		assert(m_data.at(level) != NULL);
		return *m_data.at(level);
	}

	void write(std::ostream& out) const
	{
		assert(m_data.back() != NULL);
		out << *m_data.back();
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

  protected:
	std::vector<BloomFilter*> m_data;

};

#endif
