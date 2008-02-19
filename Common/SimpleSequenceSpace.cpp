#include <algorithm>
#include "SimpleSequenceSpace.h"
#include "CommonUtils.h"

//
// Set up the 4D space to be the size of the slice passed in
//
SimpleSequenceSpace::SimpleSequenceSpace(int readLength, Coord4 minCoord, Coord4 maxCoord) : m_state(SS_LOADING)
{
	m_minCoord = minCoord;
	m_maxCoord = maxCoord;
	m_readLength = readLength;
	
	// coordinates are inclusive
	m_size.x = m_maxCoord.x - m_minCoord.x + 1;
	m_size.y = m_maxCoord.y - m_minCoord.y + 1;
	m_size.z = m_maxCoord.z - m_minCoord.z + 1;
	m_size.w = m_maxCoord.w - m_minCoord.w + 1;
	
	m_pSequences = new SequenceCollection;
}

//
// Destructor
//
SimpleSequenceSpace::~SimpleSequenceSpace()
{
	delete m_pSequences;
	m_pSequences = 0;
}

//
// Add a single read to the SimpleSequenceSpace
//
void SimpleSequenceSpace::addSequence(const PackedSeq& seq, bool boundsCheck)
{
	m_pSequences->push_back(seq);	
}

//
// Remove a read
//
void SimpleSequenceSpace::removeSequence(const PackedSeq& seq)
{
	markSequence(seq, SF_DELETE);
	
	PackedSeq& realSeq = *FindSequence(seq);
	
	// Remove this sequence as an extension to the adjacent sequences
	for(int i = 0; i <= 1; i++)
	{
		extDirection dir = (i == 0) ? SENSE : ANTISENSE;
		extDirection oppDir = oppositeDirection(dir);	
			
		for(int i = 0; i < NUM_BASES; i++)
		{	
			char currBase = BASES[i];
			// does this sequence have an extension to the deleted seq?
			if(realSeq.checkExtension(dir, currBase))
			{
				PackedSeq tempSeq(realSeq);				
				// generate the sequence that the extension is to
				char extBase = tempSeq.rotate(dir, currBase);				
				// remove the extension, this removes the reverse complement as well
				removeExtension(tempSeq, oppDir, extBase);
			}
		}
	}
}

//
// Remove the extension to this sequence from the record
//
void SimpleSequenceSpace::removeExtension(const PackedSeq& seq, extDirection dir, char base)
{
	SequenceIterPair iters = GetSequenceIterators(seq);
	removeExtensionPrivate(iters.first, dir, base);
	removeExtensionPrivate(iters.second, oppositeDirection(dir), complement(base));	
}

//
//
//
void SimpleSequenceSpace::removeExtensionPrivate(SequenceCollectionIter& seqIter, extDirection dir, char base)
{
	if(seqIter != m_pSequences->end())
	{
		seqIter->clearExtension(dir, base);	
		//seqIter->printExtension();
	}
}

//
// check if a sequence exists in the phase space
//
bool SimpleSequenceSpace::checkForSequence(const PackedSeq& seq) const
{
	assert(m_state != SS_LOADING);
	SequenceCollectionIter iter = FindSequence(seq);
	
	if(iter != m_pSequences->end())
	{
		// sequence was found
		return !iter->isFlagSet(SF_DELETE);
	}
	else
	{
		return false;
	}
}

//
//
//
void SimpleSequenceSpace::markSequence(const PackedSeq& seq, SeqFlag flag)
{
	assert(m_state == SS_READY);
	SequenceIterPair iters = GetSequenceIterators(seq);
	markSequencePrivate(iters.first, flag);
	markSequencePrivate(iters.second, flag);
}

//
//
//
void SimpleSequenceSpace::markSequencePrivate(SequenceCollectionIter& seqIter, SeqFlag flag)
{
	assert(m_state == SS_READY);
	
	if(seqIter != m_pSequences->end())
	{
		seqIter->setFlag(flag);
	}
}

//
//
//
bool SimpleSequenceSpace::checkSequenceFlag(const PackedSeq& seq, SeqFlag flag)
{
	assert(m_state == SS_READY);
	SequenceIterPair seqIters = GetSequenceIterators(seq);
	bool forwardFlag = checkSequenceFlagPrivate(seqIters.first, flag);
	bool reverseFlag = checkSequenceFlagPrivate(seqIters.second, flag);

	//assert(forwardFlag == reverseFlag);
	return (forwardFlag || reverseFlag);
}

//
//
//
bool SimpleSequenceSpace::checkSequenceFlagPrivate(SequenceCollectionIter& seqIter, SeqFlag flag)
{
	assert(m_state == SS_READY);

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
void SimpleSequenceSpace::finalizeBins()
{
	assert(m_state == SS_LOADING);
	// Sort the sequence space
	std::sort(m_pSequences->begin(), m_pSequences->end());
	
	bool checkMultiplicity = true;
	
	if(checkMultiplicity)
	{
		
		bool duplicates = checkForDuplicates();
		
		if(duplicates)
		{
			printf("duplicate sequences found, removing them\n");
			SequenceCollection temp;
			// copy the unique elements over
			std::back_insert_iterator<SequenceCollection> insertIter(temp);
			std::unique_copy(m_pSequences->begin(), m_pSequences->end(), insertIter);
			
			// swap vectors
			temp.swap(*m_pSequences);
			printf("%d sequences remain after duplicate removal\n", countAll());
		}
			
	}
	m_state = SS_FINALIZED;	
}

//
//
//
bool SimpleSequenceSpace::checkForDuplicates() const
{
	assert(m_state == SS_LOADING);
	ConstSequenceCollectionIter prev = m_pSequences->begin();
	ConstSequenceCollectionIter startIter = m_pSequences->begin() + 1;
	ConstSequenceCollectionIter endIter = m_pSequences->end();
	
	bool sorted = true;
	bool duplicates = false;
	for(ConstSequenceCollectionIter iter = startIter; iter != endIter; iter++)
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
}

//
//
//
void SimpleSequenceSpace::generateAdjacency()
{
	assert(m_state == SS_FINALIZED);
	for(SequenceCollectionIter iter = m_pSequences->begin(); iter != m_pSequences->end(); iter++)
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

//
// Calculate the extension of this sequence in the direction given
//

/* OLDE
HitRecord SimpleSequenceSpace::calculateExtension(const PackedSeq& currSeq, extDirection dir) const
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

HitRecord SimpleSequenceSpace::calculateExtension(const PackedSeq& currSeq, extDirection dir) const
{	
	
	// Create the return structure
	HitRecord hitRecord;
	
	// Get the iterators to the sequence and its reverse complement	
	SequenceIterPair iters = GetSequenceIterators(currSeq);
	
	extDirection oppDir = oppositeDirection(dir);
	
	// Check for the existance of the 4 possible extensions
	for(int i  = 0; i < NUM_BASES; i++)
	{
		char currBase = BASES[i];
		char compBase = complement(currBase);
		
		bool fwdExt = checkSequenceExtensionPrivate(iters.first, dir, currBase);
		bool revExt = checkSequenceExtensionPrivate(iters.second, oppDir, compBase);
			
		// Does this sequence have an extension?
		if(fwdExt || revExt)
		{
			PackedSeq extSeq(currSeq);
			extSeq.rotate(dir, currBase);
			
			// is there a forward extension?
			if(fwdExt)
			{
				hitRecord.addHit(extSeq, false);
			}
			else
			{
				// extension is of the reverse complement
				hitRecord.addHit(extSeq, true);	
			}
		}
	}
		
	return hitRecord;
}

//
//
//
bool SimpleSequenceSpace::hasParent(const PackedSeq& seq) const
{
	assert(m_state == SS_READY);
	SequenceIterPair iters = GetSequenceIterators(seq);
	bool forwardFlag = hasParentPrivate(iters.first);
	bool reverseFlag = hasChildPrivate(iters.second);
	
	// assert that the sequence and its reverse complement have identical flags
	//assert(forwardFlag == reverseFlag);
	return (forwardFlag || reverseFlag);
}

//
//
//
bool SimpleSequenceSpace::hasParentPrivate(SequenceCollectionIter seqIter) const
{
	assert(m_state == SS_READY);
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
bool SimpleSequenceSpace::hasChild(const PackedSeq& seq) const
{
	assert(m_state == SS_READY);
	SequenceIterPair iters = GetSequenceIterators(seq);
	bool forwardFlag = hasChildPrivate(iters.first);
	bool reverseFlag = hasParentPrivate(iters.second);

	// assert that the sequence and its reverse complement have identical flags
	//assert(forwardFlag == reverseFlag);
	return (forwardFlag || reverseFlag);
}

//
//
//
bool SimpleSequenceSpace::hasChildPrivate(SequenceCollectionIter seqIter) const
{
	assert(m_state == SS_READY);
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
bool SimpleSequenceSpace::checkSequenceExtensionPrivate(SequenceCollectionIter& seqIter, extDirection dir, char base) const
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
// get the multiplicity of the sequence
//
int SimpleSequenceSpace::getMultiplicity(const PackedSeq& seq)
{
	assert(false);
	return 0;
}

//
// print every read's multiplicity
//
void SimpleSequenceSpace::printAll() const
{
	assert(false);
}

int SimpleSequenceSpace::countAll() const
{
	return m_pSequences->size();
}

//
// Get the iterators pointing to the sequence and its reverse complement
//
SequenceIterPair SimpleSequenceSpace::GetSequenceIterators(const PackedSeq& seq) const
{
	SequenceIterPair iters;
	iters.first = FindSequence(seq);
	iters.second = FindSequence(reverseComplement(seq));
	return iters;
}

//
// Get the iterator pointing to the sequence
//
SequenceCollectionIter SimpleSequenceSpace::FindSequence(const PackedSeq& seq) const
{
	SequenceCollectionIter iter = std::lower_bound(m_pSequences->begin(), m_pSequences->end(), seq);
	if(iter != m_pSequences->end() && *iter != seq)
	{
		iter = m_pSequences->end();
	}
	return iter;
}

//
// Get the iterator pointing to the first sequence in the bin
//
SequenceCollectionIter SimpleSequenceSpace::getStartIter() const
{
	return m_pSequences->begin();
}

//
// Get the iterator pointing to the last sequence in the bin
//
SequenceCollectionIter SimpleSequenceSpace::getEndIter() const
{
	return m_pSequences->end();
}



