#include "PhaseSpace.h"


// Set up the 4D vectors to all be the readlength
PhaseSpace::PhaseSpace(int readLength) : m_phaseSpace(readLength, Bin3D(readLength, Bin2D(readLength, Bin1D(readLength))))
{
	
	
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
	m_phaseSpace[c.x][c.y][c.z][c.w][seq]++;
}

// check if a sequence exists in the phase space
bool PhaseSpace::checkForSequence(const Sequence& seq) const
{
	// calculate the coordinate for the sequence
	Coord4 c = PhaseSpace::SequenceToCoord4(seq);
	
	// check for the existance of the sequence
	if(	m_phaseSpace[c.x][c.y][c.z][c.w].count(seq) > 0)
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

bool PhaseSpace::hasParent(const Sequence& seq)
{
	HitRecord parents = calculateExtension(seq, ANTISENSE);
	return (parents.getNumHits() > 0);
}

bool PhaseSpace::hasChild(const Sequence& seq)
{
	HitRecord children = calculateExtension(seq, SENSE);
	return (children.getNumHits() > 0);
}

// get the multiplicity of the sequence
int PhaseSpace::getMultiplicity(const Sequence& seq, const Coord4& c)
{
	return m_phaseSpace[c.x][c.y][c.z][c.w][seq];
}

// print every read's multiplicity
void PhaseSpace::printAll() const
{
	// hideous nested loop
	for(Bin4D::const_iterator iter4 = m_phaseSpace.begin(); iter4 != m_phaseSpace.end(); iter4++)
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
	
	std::map<std::string, int> baseCount;
	
	std::string twomer;
	for(int i = 0; i < strLen - 1; i++)
	{
		twomer = seq.substr(i, 2);
		baseCount[twomer]++;
	}
	
	Coord4 c;
	
	c.x = abs(  2*(baseCount["CC"] + baseCount["CA"] + baseCount["AA"] + baseCount["AC"] + baseCount["TC"] + baseCount["AG"]) + (baseCount["CG"] + baseCount["GC"] + baseCount["AT"] + baseCount["TA"]) - strLen + 1);
	c.y = baseCount["CC"] + baseCount["GG"] + baseCount["AC"] + baseCount["GT"] + baseCount["CG"] + baseCount["GC"] + baseCount["AG"] + baseCount["CT"];
	c.z = baseCount["CC"] + baseCount["GG"] + baseCount["CA"] + baseCount["TG"] + baseCount["CG"] + baseCount["GC"] + baseCount["TC"] + baseCount["GA"];
	c.w = baseCount["CC"] + baseCount["GG"] + baseCount["AC"] + baseCount["GT"] + baseCount["CA"] + baseCount["TG"] + baseCount["AA"] + baseCount["TT"];
	
	return c;
}
