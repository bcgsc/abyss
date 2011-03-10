#include "config.h"
#include "SequenceCollection.h"
#include "Log.h"
#include "Common/Options.h"
#include "StringUtil.h" // for toSI
#include "Timer.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <unistd.h> // for sbrk

using namespace std;

SequenceCollectionHash::SequenceCollectionHash()
	: m_seqObserver(NULL), m_adjacencyLoaded(false)
{
#if HAVE_GOOGLE_SPARSE_HASH_MAP
	// sparse_hash_set uses 2.67 bits per element on a 64-bit
	// architecture and 2 bits per element on a 32-bit architecture.
	// The number of elements is rounded up to a power of two.
	if (opt::rank >= 0) {
		// Make room for 200 million k-mers. Approximately 58 million
		// 96-mers fit into 2 GB of ram, which results in a hash load
		// of 0.216, and approximately 116 million 32-mers, which
		// results in a hash load of 0.432.
		m_data.rehash(200000000);
		m_data.min_load_factor(0.2);
	} else {
		// Allocate a big hash for a single processor.
		m_data.rehash(1<<29);
		m_data.max_load_factor(0.4);
	}
#endif
}

/** Add the specified k-mer to this collection. */
void SequenceCollectionHash::add(const Kmer& seq)
{
	bool rc;
	SequenceCollectionHash::iterator it = find(seq, rc);
	if (it == m_data.end()) {
#if HAVE_GOOGLE_SPARSE_HASH_MAP
		if (m_data.size() == 0) {
			/* sparse_hash_set requires that set_deleted_key()
			 * is called before calling erase(). */
			Kmer rc(reverseComplement(seq));
			assert(rc != seq);
			m_data.set_deleted_key(rc);
		}
#endif
		m_data.insert(make_pair(seq, KmerData()));
	} else
		it->second.addMultiplicity(rc ? ANTISENSE : SENSE);
}

/** Clean up by erasing sequences flagged as deleted.
 * @return the number of sequences erased
 */
unsigned SequenceCollectionHash::cleanup()
{
	Timer(__func__);
	unsigned count = 0;
	for (iterator it = m_data.begin(); it != m_data.end();) {
		if (it->second.deleted()) {
			m_data.erase(it++);
			count++;
		} else
			++it;
	}
	shrink();
	return count;
}

/** Return the complement of the specified base.
 * If the assembly is in colour space, this is a no-op.
 */
static inline uint8_t complementBaseCode(uint8_t base)
{
	return opt::colourSpace ? base : ~base & 0x3;
}

/** Add an edge to this k-mer. */
bool SequenceCollectionHash::setBaseExtension(
		const Kmer& kmer, extDirection dir, uint8_t base)
{
	bool rc;
	SequenceCollectionHash::iterator it = find(kmer, rc);
	if (it == m_data.end())
		return false;
	bool palindrome = kmer.isPalindrome();
	if (!rc || palindrome)
		it->second.setBaseExtension(dir, base);
	if (rc || palindrome)
		it->second.setBaseExtension(!dir, complementBaseCode(base));
	return true;
}

/** Remove the specified extensions from this k-mer. */
void SequenceCollectionHash::removeExtension(const Kmer& kmer,
		extDirection dir, SeqExt ext)
{
	bool rc;
	SequenceCollectionHash::iterator it = find(kmer, rc);
	assert(it != m_data.end());
	bool palindrome = kmer.isPalindrome();
	if (!rc || palindrome)
		it->second.removeExtension(dir, ext);
	if (rc || palindrome)
		it->second.removeExtension(!dir, ~ext);
	notify(*it);
}

void SequenceCollectionHash::setFlag(const Kmer& key, SeqFlag flag)
{
	bool rc;
	SequenceCollectionHash::iterator it = find(key, rc);
	assert(it != m_data.end());
	it->second.setFlag(rc ? complement(flag) : flag);
}

void SequenceCollectionHash::wipeFlag(SeqFlag flag)
{
	for (SequenceCollectionHash::iterator it = m_data.begin();
			it != m_data.end(); ++it)
		it->second.clearFlag(flag);
}

/** Start of the data segment. */
static intptr_t sbrk0 = reinterpret_cast<intptr_t>(sbrk(0));

/** Print the load of the hash table. */
void SequenceCollectionHash::printLoad() const
{
	ptrdiff_t bytes = reinterpret_cast<intptr_t>(sbrk(0)) - sbrk0;
	size_t size = m_data.size();
	size_t buckets = m_data.bucket_count();
	logger(1) << "Hash load: " << size << " / " << buckets << " = "
		<< setprecision(3) << (float)size / buckets
		<< " using " << toSI(bytes) << "B" << endl;
}

/** Return an iterator pointing to the specified k-mer or its
 * reverse complement. Return in rc whether the sequence is reversed.
 */
SequenceCollectionHash::iterator SequenceCollectionHash::find(
		const Kmer& key, bool& rc)
{
	SequenceCollectionHash::iterator it = find(key);
	if (it != m_data.end()) {
		rc = false;
		return it;
	} else {
		rc = true;
		return find(reverseComplement(key));
	}
}

/** Return an iterator pointing to the specified k-mer or its
 * reverse complement. Return in rc whether the sequence is reversed.
 */
SequenceCollectionHash::const_iterator SequenceCollectionHash::find(
		const Kmer& key, bool& rc) const
{
	SequenceCollectionHash::const_iterator it = find(key);
	if (it != m_data.end()) {
		rc = false;
		return it;
	} else {
		rc = true;
		return find(reverseComplement(key));
	}
}

/** Return the sequence and data of the specified key.
 * The key sequence may not contain data. The returned sequence will
 * contain data.
 */
const SequenceCollectionHash::value_type& SequenceCollectionHash::
getSeqAndData(const Kmer& key) const
{
	bool rc;
	SequenceCollectionHash::const_iterator it = find(key, rc);
	// rc should not be ignored. This seems quite dubious.
	// The edges of this k-mer should be complemented.
	assert(it != m_data.end());
	return *it;
}

/** Return the data of the specified key. */
bool SequenceCollectionHash::getSeqData(const Kmer& key,
		ExtensionRecord& extRecord, int& multiplicity) const
{
	bool rc;
	SequenceCollectionHash::const_iterator it = find(key, rc);
	if (it == m_data.end())
		return false;
	const KmerData data = it->second;
	extRecord = rc ? ~data.extension() : data.extension();
	multiplicity = data.getMultiplicity();
	return true;
}

#include <cstdio>

/** Write this collection to disk.
 * @param path does not include the extension
 */
void SequenceCollectionHash::store(const char* path)
{
	assert(path != NULL);
#if HAVE_GOOGLE_SPARSE_HASH_MAP
	ostringstream s;
	s << path;
	if (opt::rank >= 0)
		s << '-' << setfill('0') << setw(3) << opt::rank;
	s << ".kmer";
	FILE* f = fopen(s.str().c_str(), "w");
	if (f == NULL) {
		perror(s.str().c_str());
		exit(EXIT_FAILURE);
	}
	shrink();
	m_data.write_metadata(f);
	m_data.write_nopointer_data(f);
	fclose(f);
#else
	// Not supported.
	assert(false);
	exit(EXIT_FAILURE);
#endif
}

/** Load this collection from disk. */
void SequenceCollectionHash::load(const char* path)
{
#if HAVE_GOOGLE_SPARSE_HASH_MAP
	FILE* f = fopen(path, "r");
	if (f == NULL) {
		perror(path);
		exit(EXIT_FAILURE);
	}
	m_data.read_metadata(f);
	m_data.read_nopointer_data(f);
	fclose(f);
	m_adjacencyLoaded = true;
#else
	(void)path;
	// Not supported.
	assert(false);
	exit(EXIT_FAILURE);
#endif
}

/** Indicate that this is a colour-space collection. */
void SequenceCollectionHash::setColourSpace(bool flag)
{
	if (!m_data.empty())
		assert(opt::colourSpace == flag);
	opt::colourSpace = flag;
}
