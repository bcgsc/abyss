#include <algorithm>

#if 0

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
	m_pSequences = new SequenceDataHash(10000000);
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
	m_pSequences->insert(seq);
}

//
// Remove a read
//
void SequenceCollectionHash::remove(const PackedSeq& seq)
{
	setFlag(seq, SF_DELETE);
	setFlag(reverseComplement(seq), SF_DELETE);	
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
// Remove the extension to this sequence from the record
//
void SequenceCollectionHash::removeExtension(const PackedSeq& seq, extDirection dir, char base)
{
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	removeExtensionByIter(iters.first, dir, base);
	removeExtensionByIter(iters.second, oppositeDirection(dir), complement(base));	
}

//
//
//
void SequenceCollectionHash::removeExtensionByIter(SequenceCollectionHashIter& seqIter, extDirection dir, char base)
{
	if(seqIter != m_pSequences->end())
	{
		const_cast<PackedSeq&>(*seqIter).clearExtension(dir, base);	
		//seqIter->printExtension();
	}
}

//
// check if a sequence exists in the phase space
//
ResultPair SequenceCollectionHash::exists(const PackedSeq& seq) const
{
	assert(m_state != CS_LOADING);
	ResultPair rp;
	
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	rp.forward = existsByIter(iters.first);
	rp.reverse = existsByIter(iters.second);
	return rp;
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
	/*
	assert(m_state == CS_LOADING);
	// Sort the sequence space
	std::sort(m_pSequences->begin(), m_pSequences->end());
	
	bool checkMultiplicity = true;
	
	if(checkMultiplicity)
	{
		
		bool duplicates = checkForDuplicates();
		
		if(duplicates)
		{
			printf("duplicate sequences found, removing them\n");
			SequenceData temp;
			// copy the unique elements over
			std::back_insert_iterator<SequenceData> insertIter(temp);
			std::unique_copy(m_pSequences->begin(), m_pSequences->end(), insertIter);
			
			// swap vectors
			temp.swap(*m_pSequences);
			printf("%d sequences remain after duplicate removal\n", count());
		}
			
	}
	*/
	m_state = CS_FINALIZED;
}

//
//
//
bool SequenceCollectionHash::checkForDuplicates() const
{
	/*
	assert(m_state == CS_LOADING);
	ConstSequenceCollectionHashIter prev = m_pSequences->begin();
	ConstSequenceCollectionHashIter startIter = m_pSequences->begin() + 1;
	ConstSequenceCollectionHashIter endIter = m_pSequences->end();
	
	bool sorted = true;
	bool duplicates = false;
	for(ConstSequenceCollectionHashIter iter = startIter; iter != endIter; iter++)
	{
		bool equal = false;
		if(*prev == *iter)
		{
			equal = true;
			duplicates = true;
		}
				
		if(!(*prev < *iter || equal))
		{
			sorted = false;
			printf("sort failure: %s < %s\n", prev->decode().c_str(), iter->decode().c_str());
			break;
		}
		
		if(*prev == *iter)
		{
			duplicates = true;
		}
		
		prev = iter;
	}
	
	assert(sorted);
	
	return duplicates;
	*/
	return false;
}

/*
//
//
//
void SequenceCollection::generateAdjacency()
{
	assert(m_state == SS_FINALIZED);
	for(SequenceCollectionHashIter iter = m_pSequences->begin(); iter != m_pSequences->end(); iter++)
	{
		for(int i = 0; i <= 1; i++)
		{
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			for(int j = 0; j < NUM_BASES; j++)
			{
				char currBase = BASES[j];
				PackedSeq testSeq(*iter);
				testSeq.rotate(dir, currBase);
				PackedSeq rc = reverseComplement(testSeq);
				
				if(checkForSequence(testSeq) || checkForSequence(rc))
				{
					iter->setExtension(dir, currBase);
				}
			}
		}
		
		//iter->printExtension();
	}
	m_state = SS_READY;	
	printf("done generating adjacency\n");
}
*/

//
// Calculate the extension of this sequence in the direction given
//

/* OLDE
HitRecord SequenceCollection::calculateExtension(const PackedSeq& currSeq, extDirection dir) const
{	
	PSequenceVector extVec;
	makeExtensions(currSeq, dir, extVec);

	// Create the return structure
	HitRecord hitRecord;
	// test for all the extensions of this sequence
	for(ConstPSequenceVectorIterator iter = extVec.begin(); iter != extVec.end(); iter++)
	{	
		// Todo: clean this up
		const PackedSeq& seq = *iter;
		PackedSeq rcSeq = reverseComplement(seq);
		
		if(checkForSequence(seq))
		{
			hitRecord.addHit(seq, false);
		}
		else if(checkForSequence(rcSeq))
		{
			hitRecord.addHit(seq, true);
		}	
	}
	
	return hitRecord;
}
*/

//
//
//
bool SequenceCollectionHash::hasParent(const PackedSeq& seq) const
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
bool SequenceCollectionHash::hasChild(const PackedSeq& seq) const
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
ResultPair SequenceCollectionHash::checkExtension(const PackedSeq& seq, extDirection dir, char base) const
{
	assert(m_state == CS_FINALIZED);
	ResultPair rp;
	SequenceHashIterPair iters = GetSequenceIterators(seq);
	rp.forward = checkExtensionByIter(iters.first, dir, base);
	rp.reverse = checkExtensionByIter(iters.second, oppositeDirection(dir), complement(base));
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
SequenceCollectionHashIter SequenceCollectionHash::getStartIter()
{
	return m_pSequences->begin();
}

//
// Get the iterator pointing to the last sequence in the bin
//
SequenceCollectionHashIter SequenceCollectionHash::getEndIter()
{
	return m_pSequences->end();
}

#endif
