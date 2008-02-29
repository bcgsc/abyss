#include <algorithm>
#include "AssemblyData.h"
#include "CommonUtils.h"
#if 0 
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
#endif