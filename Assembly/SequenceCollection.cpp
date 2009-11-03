#include "SequenceCollection.h"
#include "Log.h"
#include "Common/Options.h"
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
#if HAVE_GOOGLE_SPARSE_HASH_SET
	// sparse_hash_set uses 2.67 bits per element on a 64-bit
	// architecture and 2 bits per element on a 32-bit architecture.
	// The number of elements is rounded up to a power of two.
	if (opt::rank >= 0) {
		// Make room for 200 million k-mers. Approximately 58 million
		// 96-mers fit into 2 GB of ram, which results in a hash load
		// of 0.216, and approximately 116 million 32-mers, which
		// results in a hash load of 0.432.
		m_pSequences = new SequenceDataHash(200000000);
	} else {
		// Allocate a big hash for a single processor.
		m_pSequences = new SequenceDataHash(1<<29);
		m_pSequences->max_load_factor(0.4);
	}

	/* sparse_hash_set requires you call set_deleted_key() before
	* calling erase().
	* See http://google-sparsehash.googlecode.com
	* /svn/trunk/doc/sparse_hash_set.html#4
	*/
	Kmer deleted_key;
	memset(&deleted_key, 0xff, sizeof deleted_key);
	m_pSequences->set_deleted_key(deleted_key);
#else
	m_pSequences = new SequenceDataHash();
#endif
}

//
// Destructor
//
SequenceCollectionHash::~SequenceCollectionHash()
{
	delete m_pSequences;
	m_pSequences = 0;
}

/** Add the specified k-mer to this collection. */
void SequenceCollectionHash::add(const Kmer& seq)
{
	bool rc;
	SequenceCollectionHash::iterator it = find(seq, rc);
	if (it == m_pSequences->end())
		m_pSequences->insert(make_pair(seq, KmerData()));
	else
		it->second.addMultiplicity(rc ? ANTISENSE : SENSE);
}

/** Clean up by erasing sequences flagged as deleted.
 * @return the number of sequences erased
 */
unsigned SequenceCollectionHash::cleanup()
{
	Timer(__func__);
	unsigned count = 0;
	for (iterator it = m_pSequences->begin();
			it != m_pSequences->end();) {
		if (it->second.deleted()) {
			m_pSequences->erase(it++);
			count++;
		} else
			++it;
	}
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
		const Kmer& seq, extDirection dir, uint8_t base)
{
	iteratorPair iters = findBoth(seq);
	return setBaseExtensionByIter(iters.first, dir, base)
		|| setBaseExtensionByIter(iters.second,
				oppositeDirection(dir), complementBaseCode(base));
}

bool SequenceCollectionHash::setBaseExtensionByIter(
		iterator seqIter, extDirection dir, uint8_t base)
{
	if(seqIter != m_pSequences->end())
	{
		seqIter->second.setBaseExtension(dir, base);
		return true;
		//seqIter->printExtension();
	}	
	return false;
}

/** Remove the specified extensions from this k-mer. */
void SequenceCollectionHash::removeExtension(const Kmer& seq,
		extDirection dir, SeqExt ext)
{
	iteratorPair iters = findBoth(seq);
	bool found = removeExtensionByIter(iters.first, dir, ext)
		|| removeExtensionByIter(iters.second, !dir, ~ext);
	assert(found);
	(void)found;
	notify(getSeqAndData(iters));
}

bool SequenceCollectionHash::removeExtensionByIter(
		iterator seqIter, extDirection dir, SeqExt ext)
{
	if (seqIter != m_pSequences->end()) {
		seqIter->second.removeExtension(dir, ext);
		return true;
	} else
		return false;
}

void SequenceCollectionHash::setFlag(const Kmer& key,
		SeqFlag flag)
{
	bool rc;
	SequenceCollectionHash::iterator it = find(key, rc);
	assert(it != m_pSequences->end());
	it->second.setFlag(flag);
}

void SequenceCollectionHash::wipeFlag(SeqFlag flag)
{
	for (SequenceCollectionHash::iterator it = m_pSequences->begin();
			it != m_pSequences->end(); ++it)
		it->second.clearFlag(flag);
}

/** Print the load of the hash table. */
void SequenceCollectionHash::printLoad() const
{
	size_t size = m_pSequences->size();
	size_t buckets = m_pSequences->bucket_count();
	PrintDebug(1, "Hash load: %zu / %zu = %f\n",
			size, buckets, (float)size / buckets);
}

/** Return the iterators pointing to the specified k-mer and its
 * reverse complement.
 */
SequenceCollectionHash::iteratorPair SequenceCollectionHash::findBoth(
		const Kmer& seq)
{
	iteratorPair iters;
	iters.first = find(seq);
	if (iters.first != m_pSequences->end())
		iters.second = m_pSequences->end();
	else
		iters.second = find(reverseComplement(seq));
	return iters;
}

/** Return the sequence and data of the specified iterator pair. */
const SequenceCollectionHash::value_type& SequenceCollectionHash::
getSeqAndData(const iteratorPair& iters) const
{
	if (iters.first != m_pSequences->end())
		return *iters.first;
	if (iters.second != m_pSequences->end())
		return *iters.second;
	assert(false);
	exit(EXIT_FAILURE);
}

/** Return an iterator pointing to the specified k-mer. */
SequenceCollectionHash::iterator SequenceCollectionHash::find(
		const Kmer& key)
{
	return m_pSequences->find(key);
}

/** Return an iterator pointing to the specified k-mer. */
SequenceCollectionHash::const_iterator SequenceCollectionHash::find(
		const Kmer& key) const
{
	return m_pSequences->find(key);
}

/** Return an iterator pointing to the specified k-mer or its
 * reverse complement. Return in rc whether the sequence is reversed.
 */
SequenceCollectionHash::iterator SequenceCollectionHash::find(
		const Kmer& key, bool& rc)
{
	SequenceCollectionHash::iterator it = find(key);
	if (it != m_pSequences->end()) {
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
	if (it != m_pSequences->end()) {
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
	assert(it != m_pSequences->end());
	return *it;
}

/** Return the data of the specified key. */
bool SequenceCollectionHash::getSeqData(const Kmer& key,
		ExtensionRecord& extRecord, int& multiplicity) const
{
	bool rc;
	SequenceCollectionHash::const_iterator it = find(key, rc);
	if (it == m_pSequences->end())
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
void SequenceCollectionHash::store(const char* path) const
{
	assert(path != NULL);
#if HAVE_GOOGLE_SPARSE_HASH_SET
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
	m_pSequences->resize(0); // Shrink the hash table.
	m_pSequences->write_metadata(f);
	m_pSequences->write_nopointer_data(f);
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
#if HAVE_GOOGLE_SPARSE_HASH_SET
	FILE* f = fopen(path, "r");
	if (f == NULL) {
		perror(path);
		exit(EXIT_FAILURE);
	}
	m_pSequences->read_metadata(f);
	m_pSequences->read_nopointer_data(f);
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
	if (m_pSequences->size() > 0)
		assert(opt::colourSpace == flag);
	opt::colourSpace = flag;
}
