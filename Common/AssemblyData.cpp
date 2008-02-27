#include <algorithm>
#include "AssemblyData.h"
#include "CommonUtils.h"

//
// Set up the 4D space to be the size of the slice passed in
//
AssemblyData::AssemblyData() : m_state(DS_LOADING)
{
	m_pSequences = new SequenceCollection;
}

//
// Destructor
//
AssemblyData::~AssemblyData()
{
	delete m_pSequences;
	m_pSequences = 0;
}

//
// Add a single read to the AssemblyData
//
void AssemblyData::addSequence(const PackedSeq& seq)
{
	m_pSequences->add(seq);	
}

//
// Remove a read
//
void AssemblyData::removeSequence(const PackedSeq& seq)
{
	// This removes the reverse complement as well
	m_pSequences->remove(seq);
	
	// Remove this sequence as an extension to the adjacent sequences
	for(int i = 0; i <= 1; i++)
	{
		extDirection dir = (i == 0) ? SENSE : ANTISENSE;
		extDirection oppDir = oppositeDirection(dir);	
			
		for(int i = 0; i < NUM_BASES; i++)
		{	
			char currBase = BASES[i];
			// does this sequence have an extension to the deleted seq?
			ResultPair hasExt  = m_pSequences->checkExtension(seq, dir, currBase);
			if(hasExt.forward || hasExt.reverse)
			{
				PackedSeq tempSeq(seq);	
				// generate the sequence that the extension is to
				char extBase = tempSeq.rotate(dir, currBase);				
				// remove the extension, this removes the reverse complement as well
				m_pSequences->removeExtension(tempSeq, oppDir, extBase);
			}
		}
	}
}


//
// check if a sequence exists in the phase space
//
bool AssemblyData::checkForSequence(const PackedSeq& seq) const
{
	ResultPair exists = m_pSequences->exists(seq);
	return exists.forward || exists.reverse;
}

//
//
//
void AssemblyData::markSequence(const PackedSeq& seq, SeqFlag flag)
{
	m_pSequences->setFlag(seq, flag);
}

//
//
//
bool AssemblyData::checkSequenceFlag(const PackedSeq& seq, SeqFlag flag)
{
	ResultPair flagSet = m_pSequences->checkFlag(seq, flag);
	return flagSet.forward || flagSet.reverse;
}

//
//
//
void AssemblyData::finalize()
{
	m_pSequences->finalize();
}

//
//
//
void AssemblyData::generateAdjacency()
{
	assert(m_state == DS_LOADING);
	int count = 0;
	SequenceCollectionIter endIter  = m_pSequences->getEndIter();
	for(SequenceCollectionIter iter = m_pSequences->getStartIter(); iter != endIter; ++iter)
	{
		if(count % 100000 == 0)
		{
			printf("curr: %d\n", countAll());
			printf("generated for %d\n", count);
		}
		count++;
		
		for(int i = 0; i <= 1; i++)
		{
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			SeqExt extension;
			for(int j = 0; j < NUM_BASES; j++)
			{
				char currBase = BASES[j];
				PackedSeq testSeq(*iter);
				testSeq.rotate(dir, currBase);
				
				if(checkForSequence(testSeq))
				{
					extension.SetBase(currBase);
				}
			}
			m_pSequences->setExtension(*iter, dir, extension);			
		}
		
		
		//iter->printExtension();
	}
	
	m_state = DS_READY;
	printf("done generating adjacency\n");
}

//
// Calculate the extension of this sequence in the direction given
//

/* OLDE
HitRecord AssemblyData::calculateExtension(const PackedSeq& currSeq, extDirection dir) const
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

HitRecord AssemblyData::calculateExtension(const PackedSeq& currSeq, extDirection dir) const
{	
	
	// Create the return structure
	HitRecord hitRecord;
	
	// Check for the existance of the 4 possible extensions
	for(int i  = 0; i < NUM_BASES; i++)
	{
		char currBase = BASES[i];
		
		// SLOW
		ResultPair hasExt = m_pSequences->checkExtension(currSeq, dir, currBase);
			
		// Does this sequence have an extension?
		if(hasExt.forward || hasExt.reverse)
		{
			PackedSeq extSeq(currSeq);
			extSeq.rotate(dir, currBase);
			
			// is there a forward extension?
			if(hasExt.forward)
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
bool AssemblyData::hasParent(const PackedSeq& seq) const
{
	return m_pSequences->hasParent(seq);
}


//
//
//
bool AssemblyData::hasChild(const PackedSeq& seq) const
{
	return m_pSequences->hasChild(seq);
}

//
//
//
SequenceCollectionIter AssemblyData::getStartIter() const
{
	return m_pSequences->getStartIter();
}

//
//
//
SequenceCollectionIter AssemblyData::getEndIter() const
{
	return m_pSequences->getEndIter();
}

//
//
//
int AssemblyData::countAll() const
{
	return m_pSequences->count();
}
