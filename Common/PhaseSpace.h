#ifndef PHASESPACE_H
#define PHASESPACE_H

#include <vector>
#include <ext/hash_map>
#include <stdio.h>
#include "Sequence.h"
#include "PackedSeq.h"
#include "HitRecord.h"
#include "CommonDefs.h"

using namespace std;
using namespace __gnu_cxx;

typedef hash_map<Sequence, int> BinItem;
typedef std::vector<BinItem> Bin1D;
typedef std::vector<Bin1D> Bin2D;
typedef std::vector<Bin2D> Bin3D;
typedef std::vector<Bin3D> Bin4D;

enum PointClassification
{
	PC_BORDER,
	PC_INTERNAL,
	PC_INVALID
};

struct Coord4
{
	int x;
	int y;
	int z;
	int w;
};

class PhaseSpace
{
	public:
	
		//Allocates phase space
		PhaseSpace(int readLength, Coord4 startCoord, Coord4 size);
		
		//Deallocates phase space
		~PhaseSpace();
		
		// add many sequences
		void addReads(const SequenceVector& vec);
		
		// add a sequence to a coordinate
		void addSequence(const Sequence& seq, const Coord4& c);
		
		// get the multiplicity of the sequence
		int getMultiplicity(const Sequence& seq, const Coord4& c);
		
		// check if a sequence exists
		bool checkForSequence(const Sequence& seq) const;
		
		// calculate whether this sequence has an extension in the phase space
		HitRecord calculateExtension(const Sequence& currSeq, extDirection dir) const;		
		
		// does this sequence extend from a different node?
		bool hasParent(const Sequence& seq) const;

		// does this sequence have an extension?
		bool hasChild(const Sequence& seq) const;
		
		// print everything
		void printAll() const;
		
		// compute the coordinate of a sequence
		static Coord4 SequenceToCoord4(const Sequence& seq);
		
		// compute the coordinate of a sequence
		static Coord4 SequenceToCoord4(const PackedSeq* pSeq);		

	private:
	
		static inline int base2Idx(const char c); 
		
		PhaseSpace();
		
		Bin4D* m_pPhaseSpace;
		Coord4 m_start;
		Coord4 m_size;
		int m_readLength;
};

#endif
