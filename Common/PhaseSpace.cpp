#include "PhaseSpace.h"
#include "CommonUtils.h"


// Set up the 4D space to be the size of the slice passed in
PhaseSpace::PhaseSpace(int readLength, Coord4 minCoord, Coord4 maxCoord)
{
	m_minCoord = minCoord;
	m_maxCoord = maxCoord;
	m_readLength = readLength;
	
	m_size.x = m_maxCoord.x - m_minCoord.x;
	m_size.y = m_maxCoord.y - m_minCoord.y;
	m_size.z = m_maxCoord.z - m_minCoord.z;
	m_size.w = m_maxCoord.w - m_minCoord.w;
	
	// ensure the passed in size is a hypercube, ie each dimension is equal in length
	if(m_size.x == m_size.y && m_size.x == m_size.z && m_size.x == m_size.w)
	{ 
		m_pPhaseSpace = new Bin4D(m_size.x, Bin3D(m_size.x, Bin2D(m_size.x, Bin1D(m_size.x))));
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
		PackedSeq s(*iter);
		addSequence(s);
	}
}

// Add a single read to the phasespace
void PhaseSpace::addSequence(const PackedSeq& seq)
{
	Coord4 index = SequenceToTransformCoord4(seq);
	// Bounds check
	if(index.x >= 0 && index.x < m_size.x
	&& index.y >= 0 && index.y < m_size.y
	&& index.z >= 0 && index.z < m_size.z
	&& index.w >= 0 && index.w < m_size.w)
	{		
		//printf("added %s to (%d %d %d %d)\n", seq.decode().c_str(), index.x, index.y, index.z, index.w);
		(*m_pPhaseSpace)[index.x][index.y][index.z][index.w].insert(seq);
	}
	else
	{
		Coord4 realCoord = SequenceToCoord4(seq);
		printf("sequence is out of partition! (%d %d %d %d)\n", realCoord.x, realCoord.y, realCoord.z, realCoord.w);
		assert(false);
	}
}

// check if a sequence exists in the phase space
bool PhaseSpace::checkForSequence(const PackedSeq& seq) const
{
	// calculate the coordinate for the sequence
	Coord4 index = PhaseSpace::SequenceToTransformCoord4(seq);

	// check for the existance of the sequence
	if(	(*m_pPhaseSpace)[index.x][index.y][index.z][index.w].count(seq) > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
	
}

// Calculate the extension of this sequence in the direction given
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
		PackedSeq rcSeq = seq;
		rcSeq.reverseComplement();
		  
		if(checkForSequence(seq) || checkForSequence(rcSeq))
		{
			hitRecord.addHit(seq);
		}
	}
	
	return hitRecord;

}

bool PhaseSpace::hasParent(const PackedSeq& seq) const
{
	HitRecord parents = calculateExtension(seq, ANTISENSE);
	return (parents.getNumHits() > 0);
}

bool PhaseSpace::hasChild(const PackedSeq& seq) const
{
	HitRecord children = calculateExtension(seq, SENSE);
	return (children.getNumHits() > 0);
}

// get the multiplicity of the sequence
int PhaseSpace::getMultiplicity(const PackedSeq& seq)
{
	Coord4 c = PhaseSpace::SequenceToTransformCoord4(seq);
	assert(false);
	return 0;
	//return (*m_pPhaseSpace)[c.x][c.y][c.z][c.w][seq];
}

// print every read's multiplicity
void PhaseSpace::printAll() const
{
	// hideous nested loop
	for(Bin4D::const_iterator iter4 = m_pPhaseSpace->begin(); iter4 != m_pPhaseSpace->end(); iter4++)
		for(Bin3D::const_iterator iter3 = iter4->begin(); iter3 != iter4->end(); iter3++)
			for(Bin2D::const_iterator iter2 = iter3->begin(); iter2 != iter3->end(); iter2++)
				for(Bin1D::const_iterator iter1 = iter2->begin(); iter1 != iter2->end(); iter1++)
					for(BinItem::const_iterator itemIter = iter1->begin(); itemIter != iter1->end(); itemIter++)
					{
						printf("%d %s\n", itemIter->decode().c_str());	
					}
}

// Get the iterator pointing to the first sequence in the bin
PhaseSpaceBinIter PhaseSpace::getStartIter(Coord4 c) const
{
	Coord4 transC = PhaseSpace::TransformCoordinate(c);
	return (*m_pPhaseSpace)[transC.x][transC.y][transC.z][transC.w].begin();
}
		
// Get the iterator pointing to the last sequence in the bin
PhaseSpaceBinIter PhaseSpace::getEndIter(Coord4 c) const
{
	Coord4 transC = PhaseSpace::TransformCoordinate(c);
	return (*m_pPhaseSpace)[transC.x][transC.y][transC.z][transC.w].end();
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

Coord4 PhaseSpace::TransformCoordinate(Coord4& c) const
{
	Coord4 tc;
	tc.x = c.x - m_minCoord.x;
	tc.y = c.y - m_minCoord.y;
	tc.z = c.z - m_minCoord.z;
	tc.w = c.w - m_minCoord.w;
	return tc;
}

Coord4 PhaseSpace::SequenceToTransformCoord4(const PackedSeq& seq) const
{
	Coord4 c = SequenceToCoord4(seq);
	return TransformCoordinate(c);
}

// Calculate the coordinate of this sequence
// TODO: this could be optimized
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

