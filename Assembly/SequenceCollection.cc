#include "config.h"
#include "SequenceCollection.h"
#include "Log.h"
#include "Common/Options.h"
#include "Assembly/Options.h"
#include "MemoryUtil.h"
#include "StringUtil.h" // for toSI
#include "Timer.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <sstream>

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

/** sparse_hash_set requires that set_deleted_key()
 * is called before calling erase(). This key cannot
 * be an existing kmer in m_data. This function sets
 * the deleted key and should be called after all
 * data has been loaded.
 */
void SequenceCollectionHash::setDeletedKey()
{
#if HAVE_GOOGLE_SPARSE_HASH_MAP
	for (SequenceDataHash::iterator it = m_data.begin();
			it != m_data.end(); it++) {
		key_type rc(reverseComplement(it->first));
		bool isrc;
		SequenceDataHash::iterator search = find(rc, isrc);
		// If this is false, we should have a palindrome or we're
		// doing a SS assembly.
		if (isrc || search == m_data.end()) {
			m_data.set_deleted_key(rc);
			return;
		}
	}
	logger(1) << "error: unable to set deleted key.\n";
	exit(EXIT_FAILURE);
#else
	return;
#endif
}

/** Add the specified k-mer to this collection. */
void SequenceCollectionHash::add(const key_type& seq, unsigned coverage)
{
	bool rc;
	iterator it = find(seq, rc);
	if (it == m_data.end()) {
		m_data.insert(make_pair(seq, mapped_type(SENSE, coverage)));
	} else if (coverage > 0) {
		assert(!rc || !opt::ss);
		it->second.addMultiplicity(rc ? ANTISENSE : SENSE, coverage);
	}
}

/** Clean up by erasing sequences flagged as deleted.
 * @return the number of sequences erased
 */
size_t SequenceCollectionHash::cleanup()
{
	Timer(__func__);
	size_t count = 0;
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
		const key_type& kmer, extDirection dir, uint8_t base)
{
	bool rc;
	iterator it = find(kmer, rc);
	if (it == m_data.end())
		return false;
	if (opt::ss) {
		assert(!rc);
		it->second.setBaseExtension(dir, base);
	} else {
		bool palindrome = kmer.isPalindrome();
		if (!rc || palindrome)
			it->second.setBaseExtension(dir, base);
		if (rc || palindrome)
			it->second.setBaseExtension(!dir, complementBaseCode(base));
	}
	return true;
}

/** Remove the specified extensions from this k-mer. */
void SequenceCollectionHash::removeExtension(const key_type& kmer,
		extDirection dir, SeqExt ext)
{
	bool rc;
	iterator it = find(kmer, rc);
	assert(it != m_data.end());
	if (opt::ss) {
		assert(!rc);
		it->second.removeExtension(dir, ext);
	} else {
		bool palindrome = kmer.isPalindrome();
		if (!rc || palindrome)
			it->second.removeExtension(dir, ext);
		if (rc || palindrome)
			it->second.removeExtension(!dir, ~ext);
	}
	notify(*it);
}

void SequenceCollectionHash::setFlag(const key_type& key, SeqFlag flag)
{
	bool rc;
	iterator it = find(key, rc);
	assert(it != m_data.end());
	it->second.setFlag(rc ? complement(flag) : flag);
}

void SequenceCollectionHash::wipeFlag(SeqFlag flag)
{
	for (iterator it = m_data.begin();
			it != m_data.end(); ++it)
		it->second.clearFlag(flag);
}

/** Print the load of the hash table. */
void SequenceCollectionHash::printLoad() const
{
	size_t size = m_data.size();
	size_t buckets = m_data.bucket_count();
	logger(1) << "Hash load: " << size << " / " << buckets << " = "
		<< setprecision(3) << (float)size / buckets
		<< " using " << toSI(getMemoryUsage()) << "B" << endl;
}

/** Return an iterator pointing to the specified k-mer or its
 * reverse complement. Return in rc whether the sequence is reversed.
 */
SequenceCollectionHash::iterator SequenceCollectionHash::find(
		const key_type& key, bool& rc)
{
	iterator it = find(key);
	if (opt::ss || it != m_data.end()) {
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
		const key_type& key, bool& rc) const
{
	const_iterator it = find(key);
	if (opt::ss || it != m_data.end()) {
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
getSeqAndData(const key_type& key) const
{
	bool rc;
	const_iterator it = find(key, rc);
	// rc should not be ignored. This seems quite dubious.
	// The edges of this k-mer should be complemented.
	assert(it != m_data.end());
	return *it;
}

/** Return the data of the specified key. */
bool SequenceCollectionHash::getSeqData(const key_type& key,
		ExtensionRecord& extRecord, int& multiplicity) const
{
	bool rc;
	const_iterator it = find(key, rc);
	assert(!rc || !opt::ss);
	if (it == m_data.end())
		return false;
	const mapped_type data = it->second;
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
