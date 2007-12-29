#include "PhaseSpace.h"
#include "CommonUtils.h"


// Set up the 4D space to be the size of the slice passed in
PhaseSpace::PhaseSpace(int readLength, Coord4 startCoord, Coord4 size)
{
	m_start = startCoord;
	m_size = size;
	m_readLength = readLength;
	
	// ensure the passed in size is a hypercube, ie each dimension is equal in length
	if(size.x == size.y && size.x == size.z && size.x == size.w)
	{ 
		m_pPhaseSpace = new Bin4D(size.x, Bin3D(size.x, Bin2D(size.x, Bin1D(size.x))));
	}
	else
	{
		printf("partition must be equal in each dimension\n");
		assert(false);
	}
}

// Destructor, free the memory
PhaseSpace::~PhaseSpace()
{
	delete m_pPhaseSpace;
	m_pPhaseSpace = 0;
}

// Add a vector of reads to the phase space
void PhaseSpace::addReads(const SequenceVector& vec)
{
	for(SequenceVector::const_iterator iter = vec.begin(); iter != vec.end(); iter++)
	{
		Coord4 c = SequenceToCoord4(*iter);	
		addSequence(*iter, c);
	}
}
// Add a single read to the phasespace
void PhaseSpace::addSequence(const Sequence& seq, const Coord4& c)
{
	// Bounds check
	if(c.x >= m_start.x && c.x <= m_start.x + m_size.x
	&& c.y >= m_start.y && c.y <= m_start.y + m_size.y
	&& c.z >= m_start.z && c.z <= m_start.z + m_size.z
	&& c.w >= m_start.w && c.w <= m_start.w + m_size.w)
	{
		(*m_pPhaseSpace)[c.x][c.y][c.z][c.w][seq]++;
	}
	else
	{
		printf("sequence is out of partition!\n");
		assert(false);
	}
}

// check if a sequence exists in the phase space
bool PhaseSpace::checkForSequence(const Sequence& seq) const
{
	// calculate the coordinate for the sequence
	Coord4 c = PhaseSpace::SequenceToCoord4(seq);

	// check for the existance of the sequence
	if(	(*m_pPhaseSpace)[c.x][c.y][c.z][c.w].count(seq) > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
	
}

// Calculate the extension of this sequence in the direction given
HitRecord PhaseSpace::calculateExtension(const Sequence& currSeq, extDirection dir) const
{
	// Create the extensions for this read
	SequenceVector extVec;
	makeExtensions(currSeq, dir, extVec);

	// Create the return structure
	HitRecord hitRecord;
	// test for all the extensions of this sequence
	for(ConstSequenceVectorIterator iter = extVec.begin(); iter != extVec.end(); iter++)
	{	
		if(checkForSequence(*iter) || checkForSequence(reverseComplement(*iter)))
		{
			hitRecord.addHit(*iter);
		}
	}
	
	return hitRecord;
}

bool PhaseSpace::hasParent(const Sequence& seq) const
{
	HitRecord parents = calculateExtension(seq, ANTISENSE);
	return (parents.getNumHits() > 0);
}

bool PhaseSpace::hasChild(const Sequence& seq) const
{
	HitRecord children = calculateExtension(seq, SENSE);
	return (children.getNumHits() > 0);
}

// get the multiplicity of the sequence
int PhaseSpace::getMultiplicity(const Sequence& seq, const Coord4& c)
{
	return (*m_pPhaseSpace)[c.x][c.y][c.z][c.w][seq];
}

// print every read's multiplicity
void PhaseSpace::printAll() const
{
	// hideous nested loop
	for(Bin4D::const_iterator iter4 = m_pPhaseSpace->begin(); iter4 != m_pPhaseSpace->end(); iter4++)
	{
		for(Bin3D::const_iterator iter3 = iter4->begin(); iter3 != iter4->end(); iter3++)
		{
			for(Bin2D::const_iterator iter2 = iter3->begin(); iter2 != iter3->end(); iter2++)
			{	
				for(Bin1D::const_iterator iter1 = iter2->begin(); iter1 != iter2->end(); iter1++)
				{	
					for(BinItem::const_iterator itemIter = iter1->begin(); itemIter != iter1->end(); itemIter++)
					{
						printf("%d %s\n", itemIter->second, itemIter->first.c_str());	
					}
				}		
			}
		}
	} 	
}


// Calculate the coordinate of this sequence
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

// Calculate the coordinate of this sequence
// TODO: this could be optimized
Coord4 PhaseSpace::SequenceToCoord4(const PackedSeq* pSeq)
{
	const int strLen = pSeq->getSequenceLength();
	
	int vals[4][4];
	memset(vals, 0, sizeof(int) * 4 * 4);
	
	const char* curr;	
	for(int i = 0; i < strLen - 1; i++)
	{	
		// get first base index
		int idx1 = base2Idx(pSeq->getBase(i));
		int idx2 = base2Idx(pSeq->getBase(i+1));
		
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
}
