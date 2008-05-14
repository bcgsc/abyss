#include <algorithm>

#include "SequenceCollectionHash.h"
#include "CommonUtils.h"

bool PackedSeqEqual::operator()(const PackedSeq& obj1, const PackedSeq& obj2) const
{
	//printf("comp \n");
	return obj1 == obj2;
}

size_t PackedSeqHasher::operator()(const PackedSeq& myObj) const 
{
	int code = myObj.getHashCode();
  	//printf("hash: %d\n", code);
	return code;
	//return 1;
}

//
// 
//
SequenceCollectionHash::SequenceCollectionHash() : m_state(CS_LOADING)
{
	// Initially tell the has that a lot of sequences are on the way
	// This will remove many resizes (which are very slow)
	// Maybe use an even larger value?
#if HAVE_GOOGLE_SPARSE_HASH_SET
	m_pSequences = new SequenceDataHash(2 << 28);
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

//
// Add a single read to the SequenceCollection
//
void SequenceCollectionHash::add(const PackedSeq& seq)
{	
	// Check if the sequence exists or reverse complement exists
	SequenceHashIterPair iters = GetSequenceIterators(seq);

	// If it exists of the reverse complement exists, do not add
	if(iters.first != m_pSequences->end() || iters.second != m_pSequences->end())
	{
		if(iters.first != m_pSequences->end())
		{
			const_cast<PackedSeq&>(*iters.first).addMultiplicity();
		}
		else if(iters.second != m_pSequences->end())
		{
			const_cast<PackedSeq&>(*iters.second).addMultiplicity();
		}
	}
	else
	{
		m_pSequences->insert(seq);
	}
}

//
// Remove a read
//
void SequenceCollectionHash::remove(const PackedSeq& seq)
{
	// Mark the flag as deleted and remove it from the cache of branch ends (if it is there)
	setFlag(seq, SF_DELETE);
	
	// With the reverse complement as well
	setFlag(reverseComplement(seq), SF_DELETE);
}

// get the multiplicity of a sequence
int SequenceCollectionHash::getMultiplicity(const PackedSeq& seq)
{
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	int mult = 0;
	if(iters.first != m_pSequences->end())
	{
		mult += iters.first->getMultiplicity();
	}
	
	if(iters.second != m_pSequences->end())
	{
		mult += iters.second->getMultiplicity();
	}
	assert(mult != 0);
	return mult;
}

//
// add an extension to this sequence in the record
//
void SequenceCollectionHash::setExtension(const PackedSeq& seq, extDirection dir, SeqExt extension)
{
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	setExtensionByIter(iters.first, dir, extension);
	setExtensionByIter(iters.second, oppositeDirection(dir), extension.complement());
}

//
//
//
void SequenceCollectionHash::setExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, SeqExt extension)
{
	if(seqIter != m_pSequences->end())
	{
		const_cast<PackedSeq&>(*seqIter).setExtension(dir, extension);
		//seqIter->printExtension();
	}
}

//
// Set a single base extension
//
bool SequenceCollectionHash::setBaseExtension(const PackedSeq& seq, extDirection dir, char base)
{
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	bool baseSet = false;
	baseSet = baseSet || setBaseExtensionByIter(iters.first, dir, base);
	baseSet = baseSet || setBaseExtensionByIter(iters.second,
			oppositeDirection(dir), complementBaseChar(base));
	return baseSet;
}

bool SequenceCollectionHash::setBaseExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base)
{
	if(seqIter != m_pSequences->end())
	{
		const_cast<PackedSeq&>(*seqIter).setBaseExtension(dir, base);
		return true;
		//seqIter->printExtension();
	}	
	return false;
}

//
// Remove the extension to this sequence from the record
//
void SequenceCollectionHash::removeExtension(const PackedSeq& seq, extDirection dir, char base)
{
	SequenceHashIterPair iters = GetSequenceIterators(seq);

	removeExtensionByIter(iters.first, dir, base);
	removeExtensionByIter(iters.second,
			oppositeDirection(dir), complementBaseChar(base));	
}

//
//
//
void SequenceCollectionHash::removeExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base)
{
	if(seqIter != m_pSequences->end())
	{
		const_cast<PackedSeq&>(*seqIter).clearExtension(dir, base);		
	}
}


//
// Clear the extensions for this sequence
//
void SequenceCollectionHash::clearExtensions(const PackedSeq& seq, extDirection dir)
{
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	clearExtensionsByIter(iters.first, dir);
	clearExtensionsByIter(iters.second, oppositeDirection(dir));	
}

//
//
//
void SequenceCollectionHash::clearExtensionsByIter(SequenceCollectionHashIter& seqIter, extDirection dir)
{
	if(seqIter != m_pSequences->end())
	{
		const_cast<PackedSeq&>(*seqIter).clearAllExtensions(dir);	
		//seqIter->printExtension();
	}
}
//
// check if a sequence exists in the phase space
//
bool SequenceCollectionHash::exists(const PackedSeq& seq)
{
	assert(m_state != CS_LOADING);
	ResultPair rp;
	
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	rp.forward = existsByIter(iters.first);
	rp.reverse = existsByIter(iters.second);
	return (rp.forward || rp.reverse);
}

//
// Check if this sequence exists using an iterator
//
bool SequenceCollectionHash::existsByIter(SequenceCollectionHashIter& seqIter) const
{
	if(seqIter != m_pSequences->end())
	{
		// sequence was found
		return !seqIter->isFlagSet(SF_DELETE);
	}
	else
	{
		return false;
	}	
}

//
//
//
void SequenceCollectionHash::setFlag(const PackedSeq& seq, SeqFlag flag)
{
	assert(m_state == CS_FINALIZED);
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	setFlagByIter(iters.first, flag);
	setFlagByIter(iters.second, flag);
}

//
//
//
void SequenceCollectionHash::setFlagByIter(SequenceCollectionHashIter& seqIter, SeqFlag flag)
{
	assert(m_state == CS_FINALIZED);
	
	if(seqIter != m_pSequences->end())
	{
		const_cast<PackedSeq&>(*seqIter).setFlag(flag);
	}
}

//
//
//
bool SequenceCollectionHash::checkFlag(const PackedSeq& seq, SeqFlag flag)
{
	assert(m_state == CS_FINALIZED);
	ResultPair result;
	SequenceHashIterPair seqIters = GetSequenceIterators(seq);
	result.forward = checkFlagByIter(seqIters.first, flag);
	result.reverse = checkFlagByIter(seqIters.second, flag);

	//assert(forwardFlag == reverseFlag);
	return (result.forward || result.reverse);
}

//
//
//
bool SequenceCollectionHash::checkFlagByIter(SequenceCollectionHashIter& seqIter, SeqFlag flag)
{
	assert(m_state == CS_FINALIZED);

	// Check whether the sequence and its reverse complement both have the flag set/unset
	// They SHOULD be the same and the assert will guarentee this
	if(seqIter != m_pSequences->end())
	{
		return seqIter->isFlagSet(flag);
	}
	else
	{
		return false;
	}
}

//
//
//
void SequenceCollectionHash::finalize()
{
	m_state = CS_FINALIZED;
	
	size_t num_buckets = m_pSequences->bucket_count();
	size_t num_seqs = m_pSequences->size();
	printf("hash buckets: %zu sequences: %zu load factor: %f\n", num_buckets, num_seqs, (float)num_seqs/(float)num_buckets);
	
}

//
//
//
bool SequenceCollectionHash::checkForDuplicates() const
{
	return false;
}

//
//
//
bool SequenceCollectionHash::hasParent(const PackedSeq& seq)
{
	assert(m_state == CS_FINALIZED);
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	bool forwardFlag = hasParentByIter(iters.first);
	bool reverseFlag = hasChildByIter(iters.second);
	
	// assert that the sequence and its reverse complement have identical flags
	//assert(forwardFlag == reverseFlag);
	return (forwardFlag || reverseFlag);
}

//
//
//
bool SequenceCollectionHash::hasParentByIter(SequenceCollectionHashIter seqIter) const
{
	assert(m_state == CS_FINALIZED);
	if(seqIter != m_pSequences->end())
	{
		return seqIter->hasExtension(ANTISENSE);
	}
	else
	{
		return false;
	}
}

//
//
//
bool SequenceCollectionHash::hasChild(const PackedSeq& seq)
{
	assert(m_state == CS_FINALIZED);
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	bool forwardFlag = hasChildByIter(iters.first);
	bool reverseFlag = hasParentByIter(iters.second);

	// assert that the sequence and its reverse complement have identical flags
	//assert(forwardFlag == reverseFlag);
	return (forwardFlag || reverseFlag);
}

//
//
//
bool SequenceCollectionHash::hasChildByIter(SequenceCollectionHashIter seqIter) const
{
	assert(m_state == CS_FINALIZED);
	if(seqIter != m_pSequences->end())
	{
		return seqIter->hasExtension(SENSE);
	}
	else
	{
		return false;
	}
}

//
//
//
ResultPair SequenceCollectionHash::checkExtension(const PackedSeq& seq, extDirection dir, char base)
{
	assert(m_state == CS_FINALIZED);
	ResultPair rp;
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	rp.forward = checkExtensionByIter(iters.first, dir, base);
	rp.reverse = checkExtensionByIter(iters.second,
			oppositeDirection(dir), complementBaseChar(base));
	return rp;
}

//
//
//
bool SequenceCollectionHash::checkExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base) const
{
	if(seqIter != m_pSequences->end())
	{
		return seqIter->checkExtension(dir, base);
	}
	else
	{
		return false;
	}
}

//
// get the extensions and multiplicity information for a sequence
//
bool SequenceCollectionHash::getSeqData(const PackedSeq& seq, ExtensionRecord& extRecord, int& multiplicity)
{
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	bool found = false;
	if(iters.first != m_pSequences->end())
	{
		// seq found
		extRecord.dir[SENSE] = getExtensionByIter(iters.first, SENSE);
		extRecord.dir[ANTISENSE] = getExtensionByIter(iters.first, ANTISENSE);
		multiplicity = iters.first->getMultiplicity();
		found = true;
	}
	else if(iters.second != m_pSequences->end())
	{
		// reverse seq found, swap the order of the extensions and complement them (to make sure they are in the direction of the requested seq)
		extRecord.dir[SENSE] = getExtensionByIter(iters.second, ANTISENSE).complement();
		extRecord.dir[ANTISENSE] = getExtensionByIter(iters.second, SENSE).complement();
		multiplicity = iters.second->getMultiplicity();
		found = true;
	}

	return found;
}

//
//
//
SeqExt SequenceCollectionHash::getExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir) const
{
	SeqExt ret;
	if(seqIter != m_pSequences->end())
	{
		// Sequence found
		ret = seqIter->getExtension(dir);
	}
	
	return ret;	
}

//
// Return the number of sequences held in the collection
// Note: some sequences will be marked as DELETE and will still be counted
//
int SequenceCollectionHash::count() const
{
	return m_pSequences->size();
}

//
// Get the iterators pointing to the sequence and its reverse complement
//
SequenceHashIterPair SequenceCollectionHash::GetSequenceIterators(const PackedSeq& seq) const
{
	SequenceHashIterPair iters;
	iters.first = FindSequence(seq);
	iters.second = FindSequence(reverseComplement(seq));
	return iters;
}

//
// Get the iterator pointing to the sequence
//
SequenceCollectionHashIter SequenceCollectionHash::FindSequence(const PackedSeq& seq) const
{
	return m_pSequences->find(seq);
}


//
// Get the iterator pointing to the first sequence in the bin
//
SequenceCollectionHashIter SequenceCollectionHash::getStartIter() const
{
	return m_pSequences->begin();
}

//
// Get the iterator pointing to the last sequence in the bin
//
SequenceCollectionHashIter SequenceCollectionHash::getEndIter() const
{
	return m_pSequences->end();
}

