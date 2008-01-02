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

typedef map<PackedSeq, int> BinItem;
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
		
		// add a single sequence
		void addSequence(const PackedSeq& seq);
		
		// get the multiplicity of the sequence
		int getMultiplicity(const PackedSeq& seq);
		
		// check if a sequence exists
		bool checkForSequence(const PackedSeq& seq) const;
		
		// calculate whether this sequence has an extension in the phase space
		HitRecord calculateExtension(const PackedSeq& currSeq, extDirection dir) const;		
		
		// does this sequence extend from a different node?
		bool hasParent(const PackedSeq& seq) const;

		// does this sequence have an extension?
		bool hasChild(const PackedSeq& seq) const;
		
		// print everything
		void printAll() const;
		
		// Compute the transformed (phase space) coordinate of this sequence
		Coord4 SequenceToTransformCoord4(const PackedSeq& seq) const;
		
		// compute the coordinate of a sequence
		static Coord4 SequenceToCoord4(const Sequence& seq);
		
		// compute the coordinate of a sequence
		static Coord4 SequenceToCoord4(const PackedSeq& pSeq);

	private:
	
		static inline int base2Idx(const char c); 
		
		PhaseSpace();
		
		Bin4D* m_pPhaseSpace;
		Coord4 m_start;
		Coord4 m_size;
		int m_readLength;
};

#endif
