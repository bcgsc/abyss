#ifndef COUNTINGBLOOM_H
#define COUNTINGBLOOM_H

#include <cstring>
#include <limits>
#include <cmath>
#include <vector>
#include <cassert>
#include <fstream>

using namespace std;

// Forward declaraions.
template <typename T>
class CountingBloomFilter;

// Method declarations.
template <typename T>
std::ostream& operator<<(std::ostream&, const CountingBloomFilter<T>&);

template <typename T>
class CountingBloomFilter {
public:
	CountingBloomFilter(): m_filter(NULL), m_size(0), m_sizeInBytes(0),
			    m_hashNum(0), m_kmerSize(0), m_nEntry(0),
			    m_tEntry(0), m_dFPR(0), m_countThreshold(0) { }
	CountingBloomFilter(size_t sz, unsigned hashNum, unsigned kmerSize,
			 unsigned countThreshold)
		: m_filter(new T[sz]), m_size(sz),
		  m_sizeInBytes(sz * sizeof(T)), m_hashNum(hashNum),
		  m_kmerSize(kmerSize), m_nEntry(0), m_tEntry(0), m_dFPR(0),
		  m_countThreshold(countThreshold) {
		std::memset(m_filter, 0, m_sizeInBytes);
	}
	CountingBloomFilter(const string &path);
	~CountingBloomFilter() {
		delete[] m_filter;
	}
	T operator[](size_t i) {
		return m_filter[i];
	}
	template <typename U> T minCount(const U &hashes) const {
		T min = m_filter[hashes[0] % m_size];
		for (size_t i = 1; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_size;
			if (m_filter[pos] < min)
				min = m_filter[pos];
		}
		return min;
	}
	template <typename U> bool   contains(const U &hashes) const;
	template <typename U> void   insert(const U &hashes);
	template <typename U> bool   insertAndCheck(const U &hashes);
	template <typename U> void   incrementMin(const U &hashes);
	template <typename U> void   incrementAll(const U &hashes);
	unsigned getKmerSize(void) const { return m_kmerSize; };
	unsigned getHashNum(void)  const { return m_hashNum; };
	size_t   size(void)        const { return m_size; };
	size_t   sizeInBytes(void) const { return m_sizeInBytes; };
	unsigned threshold(void)   const { return m_countThreshold; };
	size_t   popCount() const;
	size_t   popCount_threshold() const;
	double   FPR(void)  const;
        double   FPR_post_threshold(void) const;

	// Serialization interface
	struct FileHeader {
		char magic[8];
		uint32_t hlen;
		uint64_t size;
		uint32_t nhash;
		uint32_t kmer;
		double   dFPR;
		uint64_t nEntry;
		uint64_t tEntry;
	};
	void readHeader(FILE *file);
	void readFilter(const string &path);
	void writeHeader(std::ostream& out) const;
	void writeFilter(string const &path) const;
	friend std::ostream& operator<< <> (std::ostream&,
					    const CountingBloomFilter&);

private:
	// m_filter         : An array of elements of type T; the bit-array or
	//                    filter.
	// m_size           : Size of bloom filter (size of m_filter array).
	// m_sizeInBytes    : Size of the bloom filter in bytes, that is,
	//                    (m_size * sizeof(T)).
	// m_hashNum        : Number of hash functions.
	// m_kmerSize       : Size of a k-mer.
	// m_nEntry         : Number of items the bloom filter holds.
	// m_tEntry         : ?
	// m_dFPR           : Why d?
	// m_countThreshold : A count greater or equal to this threshold
	//                    establishes existence of an element in the filter.

	T        *m_filter;
	size_t   m_size;
	size_t   m_sizeInBytes;
	unsigned m_hashNum;
	unsigned m_kmerSize;
	size_t   m_nEntry;
	size_t   m_tEntry;
	double   m_dFPR;
	unsigned m_countThreshold;
};

// Method definitions

// Bloom filter operations

/*
  Use of atomic increments in incrementMin() and incrementAll():

  A atomic compare-and-swap (CAS) operation increments m_filter[pos]. The CAS
  operation takes a memory location and a value that the caller believes the
  location currently stores. If the memory location still holds that value when
  the atomic compare-and-swap executes, then a new value is stored and 'true'
  returned; otherwise, memory is left unchanged and 'false' returned.

  The value of m_filter[pos] may be changed by another thread between a read from
  that memory location and a write to it. The CAS operation is called in a loop
  until it succeeds, which ensures that a write does not happen if some other
  thread has incremented the value between this thread's read and write.

  Note that CAS operations suffer from the ABA problem.
*/


// Of the m_hashNum counters, increment all the minimum values.
template <typename T>
template <typename U>
inline void CountingBloomFilter<T>::incrementMin(const U &hashes) {
        // update flag to track if increment is done on at least one counter 
        bool updateDone=false;
	T newVal;
	T minVal = minCount(hashes);
        while (!updateDone){
            // Simple check to deal with overflow
            newVal = minVal + 1;
            if (minVal > newVal){
                return;
            }
            for (size_t i = 0; i < m_hashNum; ++i) {
                if (__sync_bool_compare_and_swap(&m_filter[hashes[i] % m_size], minVal, newVal)){
                    updateDone = true;                    
                }
            }
            // Recalculate minval because if increment fails, it needs a new minval to use and 
            // if it doesnt hava a new one, the while loop runs forever.
            if (!updateDone){
                minVal = minCount(hashes);    
            }
        }
        return;
}



// Not Used, deprecated
// Increment all the m_hashNum counters.
template <typename T>
template <typename U>
inline void CountingBloomFilter<T>::incrementAll(const U &hashes) {
	T currentVal, newVal;
	for (size_t i = 0; i < m_hashNum; ++i) {
		size_t pos = hashes[i] % m_size;
		do {
			currentVal = m_filter[pos];
			newVal     = currentVal + 1;
			if (newVal < currentVal)
				break;
		} while(!__sync_bool_compare_and_swap(&m_filter[pos],
						      currentVal, newVal));
	}
}


// Check if an element exists. If the minimum count at the m_hashNum positions
// of m_filter is more than or equal to a predefined count threshold, then the
// element is said to be present in the Bloom filter. count() therefore returns
// true when this condition is satisfied, or else, false.

template <typename T>
template <typename U>
inline bool CountingBloomFilter<T>::contains(const U &hashes) const {
	return minCount(hashes) >= m_countThreshold;
}

template <typename T>
template <typename U>
inline void CountingBloomFilter<T>::insert(const U &hashes) {
	incrementMin(hashes);
}

template <typename T>
template <typename U>
inline bool CountingBloomFilter<T>::insertAndCheck(const U &hashes) {
	bool found = contains(hashes);
	incrementMin(hashes);
	return found;
}

/* Count the number of non-zero counters. */
template <typename T>
size_t CountingBloomFilter<T>::popCount() const {
	size_t count = 0;
	for (size_t i = 0; i < m_size; ++i) {
            if (m_filter[i] != 0)
                ++count;
	}
	return count;
}


/* Count the number of above threshold counters. */
template <typename T>
size_t CountingBloomFilter<T>::popCount_threshold() const {
	size_t count = 0;
	for (size_t i = 0; i < m_size; ++i) {
	    if (m_filter[i] >= m_countThreshold )
	    ++count;
	}
	return count;
}



template <typename T>
double CountingBloomFilter<T>::FPR(void) const {
	return std::pow((double)popCount() / (double)m_size, m_hashNum);
}
template <typename T>
double CountingBloomFilter<T>::FPR_post_threshold(void) const {
	return std::pow((double)popCount_threshold() / (double)m_size, m_hashNum);
}



// Serialization interface.
template <typename T>
CountingBloomFilter<T>::CountingBloomFilter(const string &path) {
	if (m_filter != NULL)
		delete[] m_filter;
	m_filter = new T[m_size];
	readFilter(path);
}
template <typename T>
void CountingBloomFilter<T>::readFilter(const string &path) {
	FILE *fp;
	if ((fp = fopen(path.c_str(), "rb")) == NULL) {
		cerr << "ERROR: Failed to open file: " << path << "\n";
		exit(1);
	}
	readHeader(fp);
	long int lCurPos = ftell(fp);
	fseek(fp, 0, 2);
	size_t arraySizeOnDisk = ftell(fp) - sizeof(struct FileHeader);
	fseek(fp, lCurPos, 0);
	if (arraySizeOnDisk != m_sizeInBytes) {
		cerr << "ERROR: File size of " << path << " ("
		     << arraySizeOnDisk << " bytes), "
		     << "does not match size read from its header ("
		     << m_sizeInBytes << " bytes).\n";
		exit(1);
	}

	size_t nread = fread(m_filter, arraySizeOnDisk, 1, fp);
	if (nread != 1 && fclose(fp) != 0) {
		cerr << "ERROR: The bit array could not be read from the file: "
		     << path << "\n";
		exit(1);
	}
}
template <typename T>
void CountingBloomFilter<T>::readHeader(FILE *fp) {
	FileHeader header;
	char magic[9];
	if (fread(&header, sizeof(struct FileHeader), 1, fp) != 1) {
		cerr << "Failed to read header\n";
		exit(1);
	}
	memcpy(magic, header.magic, 8);
	magic[8] = '\0';
	m_size        = header.size;
	m_hashNum     = header.nhash;
	m_kmerSize    = header.kmer;
	m_sizeInBytes = m_size * sizeof(T);
}
template <typename T>
void CountingBloomFilter<T>::writeFilter(string const &path) const {
	ofstream ofs(path.c_str(), ios::out | ios::binary);
	cerr << "Writing a " << m_sizeInBytes << " byte filter to a file on disk.\n";
	ofs << *this;
	ofs.close();
	assert(ofs);
}
template <typename T>
void CountingBloomFilter<T>::writeHeader(std::ostream &out) const {
	FileHeader header;
	char magic[9];
	memcpy(header.magic, "BlOOMFXX", 8);
	memcpy(magic, header.magic, 8);
	magic[8] = '\0';
	header.hlen   = sizeof(struct FileHeader);
	header.size   = m_size;
	header.nhash  = m_hashNum;
	header.kmer   = m_kmerSize;
	header.dFPR   = m_dFPR;
	header.nEntry = m_nEntry;
	header.tEntry = m_tEntry;
	out.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));
	assert(out);
}
// Serialize the bloom filter to a C++ stream */
template <typename T>
std::ostream& operator<<(std::ostream &os, const CountingBloomFilter<T>& cbf) {
	assert(os);
	cbf.writeHeader(os);
	assert(os);
	os.write(reinterpret_cast<char*>(cbf.m_filter), cbf.m_sizeInBytes);
	assert(os);
	return os;
}

#endif // COUNTINGBLOOM_H
