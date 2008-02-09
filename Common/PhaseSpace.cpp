#include <algorithm>
#include "PhaseSpace.h"
#include "CommonUtils.h"

//
// Set up the 4D space to be the size of the slice passed in
//
PhaseSpace::PhaseSpace(int readLength, Coord4 minCoord, Coord4 maxCoord) : m_writeEnabled(true)
{
	m_minCoord = minCoord;
	m_maxCoord = maxCoord;
	m_readLength = readLength;
	
	// coordinates are inclusive
	m_size.x = m_maxCoord.x - m_minCoord.x + 1;
	m_size.y = m_maxCoord.y - m_minCoord.y + 1;
	m_size.z = m_maxCoord.z - m_minCoord.z + 1;
	m_size.w = m_maxCoord.w - m_minCoord.w + 1;
	
	m_pPhaseSpace = new Bin4D(m_size.x, Bin3D(m_size.y, Bin2D(m_size.z, Bin1D(m_size.w))));

}

//
// Destructor
//
PhaseSpace::~PhaseSpace()
{
	delete m_pPhaseSpace;
	m_pPhaseSpace = 0;
}

//
// Add a vector of reads to the phase space
//
void PhaseSpace::addReads(const SequenceVector& vec)
{
	for(SequenceVector::const_iterator iter = vec.begin(); iter != vec.end(); iter++)
	{	
		PackedSeq s(*iter);
		addSequence(s);
	}
}

//
// Add a single read to the phasespace
//
void PhaseSpace::addSequence(const PackedSeq& seq, bool boundsCheck)
{
	// Make sure the phase space is writable (ie it 
	assert(m_writeEnabled);
	
	Coord4 index = SequenceToIndex(seq);
	// Bounds check
	if(CheckValidIndex(index))
	{		
		//printf("added %s to (%d %d %d %d)\n", seq.decode().c_str(), index.x, index.y, index.z, index.w);
		
		// perform a sorted insert of the sequence into the vector
		BinItem& currBin = (*m_pPhaseSpace)[index.x][index.y][index.z][index.w];
		currBin.insert(seq);
	}
	else if(boundsCheck)
	{
		Coord4 realCoord = SequenceToCoord4(seq);
		printf("sequence is out of partition! (%d %d %d %d)\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
		assert(false);
	}
}

//
// Remove a read
//
void PhaseSpace::removeSequence(const PackedSeq& seq)
{
	Coord4 index = SequenceToIndex(seq);
	// Bounds check
	if(CheckValidIndex(index))
	{
		BinItem& currBin = (*m_pPhaseSpace)[index.x][index.y][index.z][index.w];
		currBin.erase(seq);
	}
}

//
// check if a sequence exists in the phase space
//
bool PhaseSpace::checkForSequence(const PackedSeq& seq) const
{
	// calculate the coordinate for the sequence
	Coord4 index = SequenceToIndex(seq);
	
	// Make sure we aren't requesting an invalid coordinate
	// This would indicate that the correct phase space wasn't loaded
	
	if(CheckValidIndex(index))
	{
		// Reference to the correct vector
		BinItem& currBin = (*m_pPhaseSpace)[index.x][index.y][index.z][index.w];
		// Search the SORTED vector
		return currBin.find(seq) != currBin.end();
	}
	else
	{
		Coord4 realCoord = SequenceToCoord4(seq);
		printf("sequence is out of partition! (%d %d %d %d)\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
		assert(false);		
	}
}

//
// Searches the phase space for a particular sequence and returns the reference IN the phase space to it
// this allows us to manipulate the sequences that make up the phase space (in particular mark them for deletion, etc)
//
void PhaseSpace::markSequence(const PackedSeq& seq, SeqFlag flag)
{
	// calculate the coordinate for the sequence
	Coord4 index = SequenceToIndex(seq);
	if(CheckValidIndex(index))
	{
		// Reference to the correct vector
		BinItem& currBin = (*m_pPhaseSpace)[index.x][index.y][index.z][index.w];

		PhaseSpaceBinIter itemIter = currBin.find(seq);
		if(itemIter != currBin.end())
		{
			
			assert(*itemIter == seq);
		
			// this is only valid because setting the flag doesnt not effect the relative ordering of the trees
			// const_cast is generally pretty hacky
			const_cast<PackedSeq&>(*itemIter).setFlag(flag);
		}
	}
	else
	{
		Coord4 realCoord = SequenceToCoord4(seq);
		printf("sequence is out of partition! (%d %d %d %d)\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
		assert(false);		
	}			
}

//
//
//
bool PhaseSpace::checkSequenceFlag(const PackedSeq& seq, SeqFlag flag)
{
	// calculate the coordinate for the sequence
	Coord4 index = SequenceToIndex(seq);
	if(CheckValidIndex(index))
	{
		// Reference to the correct vector
		BinItem& currBin = (*m_pPhaseSpace)[index.x][index.y][index.z][index.w];
		
		PhaseSpaceBinIter itemIter = currBin.find(seq);
		if(itemIter != currBin.end())
		{			
			return itemIter->isFlagSet(flag);
		}
		else
		{
			return false;
		}
	}
	else
	{
		Coord4 realCoord = SequenceToCoord4(seq);
		printf("sequence is out of partition! (%d %d %d %d)\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
		assert(false);		
	}
}

//
//
//
void PhaseSpace::finalizeBins(Coord4 start, Coord4 end)
{
	assert(false);
#if 0
	// Disable writes to the phase space, trim the vectors and sort them
	//m_writeEnabled = false;
	//printf("finalizing....");
	
	for(int x = start.x; x <= end.x; x++)
		for(int y = start.y; y <= end.y; y++)
			for(int z = start.z; z <= end.z; z++)
				for(int w = start.w; w <= end.w; w++)
				{
					Coord4 c = {x,y,z,w};
					Coord4 index = CoordToIndex(c);
					assert(CheckValidIndex(index));
					
				
					// Get the current bin
					BinItem& currBin = (*m_pPhaseSpace)[index.x][index.y][index.z][index.w];

					// Sort the current bin
					std::sort(currBin.begin(), currBin.end());
					
					// swap trick to downsize
					
					// what happens here is you create a new exact-sized vector which is a copy of vector1. Then swap it back into vector1
					// This allows you to trim the vector down to the minimum required size since you can't downsize an stl vector
					BinItem(currBin.begin(), currBin.end()).swap(currBin);
					//printf("after swap: %d (%d)\n", iter1->size(), iter1->capacity());
				}
#endif
}

//
// Calculate the extension of this sequence in the direction given
//
HitRecord PhaseSpace::calculateExtension(const PackedSeq& currSeq, extDirection dir) const
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

//
//
//
bool PhaseSpace::hasParent(const PackedSeq& seq) const
{
	HitRecord parents = calculateExtension(seq, ANTISENSE);
	return (parents.getNumHits() > 0);
}

//
//
//
bool PhaseSpace::hasChild(const PackedSeq& seq) const
{
	HitRecord children = calculateExtension(seq, SENSE);
	return (children.getNumHits() > 0);
}

//
// get the multiplicity of the sequence
//
int PhaseSpace::getMultiplicity(const PackedSeq& seq)
{
	assert(false);
	return 0;
}

//
// print every read's multiplicity
//
void PhaseSpace::printAll() const
{
	assert(false);
	// hideous nested loop
	for(Bin4D::const_iterator iter4 = m_pPhaseSpace->begin(); iter4 != m_pPhaseSpace->end(); iter4++)
		for(Bin3D::const_iterator iter3 = iter4->begin(); iter3 != iter4->end(); iter3++)
			for(Bin2D::const_iterator iter2 = iter3->begin(); iter2 != iter3->end(); iter2++)
				for(Bin1D::const_iterator iter1 = iter2->begin(); iter1 != iter2->end(); iter1++)
					for(BinItem::const_iterator itemIter = iter1->begin(); itemIter != iter1->end(); itemIter++)
					{
						//printf("%d %s\n", itemIter->decode().c_str());	
					}
}

int PhaseSpace::countAll() const
{
	int total = 0;
	for(Bin4D::const_iterator iter4 = m_pPhaseSpace->begin(); iter4 != m_pPhaseSpace->end(); iter4++)
		for(Bin3D::const_iterator iter3 = iter4->begin(); iter3 != iter4->end(); iter3++)
			for(Bin2D::const_iterator iter2 = iter3->begin(); iter2 != iter3->end(); iter2++)
				for(Bin1D::const_iterator iter1 = iter2->begin(); iter1 != iter2->end(); iter1++)
					total += iter1->size();
	return total;
}

//
// Get the iterator pointing to the first sequence in the bin
//
PhaseSpaceBinIter PhaseSpace::getStartIter(Coord4 c) const
{
	Coord4 index = CoordToIndex(c);
	assert(CheckValidIndex(index));
	return (*m_pPhaseSpace)[index.x][index.y][index.z][index.w].begin();
}

//
// Get the iterator pointing to the last sequence in the bin
//
PhaseSpaceBinIter PhaseSpace::getEndIter(Coord4 c) const
{
	Coord4 index = CoordToIndex(c);
	assert(CheckValidIndex(index));	
	return (*m_pPhaseSpace)[index.x][index.y][index.z][index.w].end();
}

//
// Calculate the coordinate of this sequence
//
Coord4 PhaseSpace::SequenceToCoord4(const Sequence& seq)
{
	const int strLen = seq.length();
	const char* data = seq.data();
	
	int vals[4][4];
	memset(vals, 0, sizeof(int) * 4 * 4);
	
	const char* curr;	
	for(int i = 0; i < strLen - 1; i++)
	{
		curr = data + i;
		
		// get first base index
		int idx1 = base2Idx(*curr);
		int idx2 = base2Idx(*(curr + 1));
		
		vals[idx1][idx2]++;
	}
	
	Coord4 c;
	
	int idxA = base2Idx('A');
	int idxC = base2Idx('C');
	int idxG = base2Idx('G');
	int idxT = base2Idx('T');
	
	c.x = abs(  2*(vals[idxC][idxC] + vals[idxC][idxA] + vals[idxA][idxA] + vals[idxA][idxC] + vals[idxT][idxC] + vals[idxA][idxG]) + (vals[idxC][idxG] + vals[idxG][idxC] + vals[idxA][idxT] + vals[idxT][idxA]) - strLen + 1);
	c.y = vals[idxC][idxC] + vals[idxG][idxG] + vals[idxA][idxC] + vals[idxG][idxT] + vals[idxC][idxG] + vals[idxG][idxC] + vals[idxA][idxG] + vals[idxC][idxT];
	c.z = vals[idxC][idxC] + vals[idxG][idxG] + vals[idxC][idxA] + vals[idxT][idxG] + vals[idxC][idxG] + vals[idxG][idxC] + vals[idxT][idxC] + vals[idxG][idxA];
	c.w = vals[idxC][idxC] + vals[idxG][idxG] + vals[idxA][idxC] + vals[idxG][idxT] + vals[idxC][idxA] + vals[idxT][idxG] + vals[idxA][idxA] + vals[idxT][idxT];
	
	return c;
}

//
//
//
Coord4 PhaseSpace::CoordToIndex(const Coord4& c) const
{
	Coord4 tc;
	tc.x = c.x - m_minCoord.x;
	tc.y = c.y - m_minCoord.y;
	tc.z = c.z - m_minCoord.z;
	tc.w = c.w - m_minCoord.w;
	return tc;
}

//
//
//
bool PhaseSpace::CheckValidCoordinate(const Coord4& c) const
{
	return CheckValidIndex(CoordToIndex(c));
}

//
//
//
bool PhaseSpace::CheckValidIndex(const Coord4& c) const
{
	return (c.x >= 0 && c.x < m_size.x && c.y >= 0 && c.y < m_size.y && c.z >= 0 && c.z < m_size.z && c.w >= 0 && c.w < m_size.w);
}

//
//
//
Coord4 PhaseSpace::SequenceToIndex(const PackedSeq& seq) const
{
	Coord4 c = SequenceToCoord4(seq);
	return CoordToIndex(c);
}

//
// Calculate the coordinate of this sequence
// TODO: this could be optimized
//
Coord4 PhaseSpace::SequenceToCoord4(const PackedSeq& pSeq)
{
	const int strLen = pSeq.getSequenceLength();
	
	int vals[4][4];
	memset(vals, 0, sizeof(int) * 4 * 4);
	
	for(int i = 0; i < strLen - 1; i++)
	{	
		// get first base index
		int idx1 = base2Idx(pSeq.getBase(i));
		int idx2 = base2Idx(pSeq.getBase(i+1));
		
		vals[idx1][idx2]++;
	}
	
	Coord4 c;
	
	int idxA = base2Idx('A');
	int idxC = base2Idx('C');
	int idxG = base2Idx('G');
	int idxT = base2Idx('T');
	
	c.x = abs(  2*(vals[idxC][idxC] + vals[idxC][idxA] + vals[idxA][idxA] + vals[idxA][idxC] + vals[idxT][idxC] + vals[idxA][idxG]) + (vals[idxC][idxG] + vals[idxG][idxC] + vals[idxA][idxT] + vals[idxT][idxA]) - strLen + 1);
	c.y = vals[idxC][idxC] + vals[idxG][idxG] + vals[idxA][idxC] + vals[idxG][idxT] + vals[idxC][idxG] + vals[idxG][idxC] + vals[idxA][idxG] + vals[idxC][idxT];
	c.z = vals[idxC][idxC] + vals[idxG][idxG] + vals[idxC][idxA] + vals[idxT][idxG] + vals[idxC][idxG] + vals[idxG][idxC] + vals[idxT][idxC] + vals[idxG][idxA];
	c.w = vals[idxC][idxC] + vals[idxG][idxG] + vals[idxA][idxC] + vals[idxG][idxT] + vals[idxC][idxA] + vals[idxT][idxG] + vals[idxA][idxA] + vals[idxT][idxT];
	
	return c;
}

//
//
//
int PhaseSpace::base2Idx(const char c)
{
	if(c == 'A')
	{
		return 0;
	}
	else if(c == 'C')
	{
		return 1;
	}
	else if(c == 'G')
	{
		return 2;
	}
	else if(c == 'T')
	{
		return 3;
	}
	assert(false);
	return -1;
}

