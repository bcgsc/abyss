#ifndef PHASESPACE_H
#define PHASESPACE_H

#include <vector>
#include <set>
#include <stdio.h>
#include "Sequence.h"
#include "PackedSeq.h"
#include "HitRecord.h"
#include "CommonDefs.h"

using namespace std;


typedef std::vector<PackedSeq> BinItem;
typedef std::vector<BinItem> Bin1D;
typedef std::vector<Bin1D> Bin2D;
typedef std::vector<Bin2D> Bin3D;
typedef std::vector<Bin3D> Bin4D;
typedef BinItem::iterator PhaseSpaceBinIter;
typedef BinItem::const_iterator ConstPhaseSpaceBinIter;

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
		void addSequence(const PackedSeq& seq, bool boundsCheck = false);
		
		// trim and sort the vectors
		void finalizeBins(Coord4 start, Coord4 end);
		
		// get the multiplicity of the sequence
		int getMultiplicity(const PackedSeq& seq);
		
		// check if a sequence exists
		bool checkForSequence(const PackedSeq& seq) const;
		
		// check if the coordinate is valid
		inline bool CheckValidCoordinate(const Coord4& c) const;
		
		// check if the index is valid
		inline bool CheckValidIndex(const Coord4& c) const;
		
		// calculate whether this sequence has an extension in the phase space
		HitRecord calculateExtension(const PackedSeq& currSeq, extDirection dir) const;
		
		// Get the iterator pointing to the first sequence in the bin
		PhaseSpaceBinIter getStartIter(Coord4 c) const;
		
		// Get the iterator pointing to the last sequence in the bin
		PhaseSpaceBinIter getEndIter(Coord4 c) const;
		
		// does this sequence extend from a different node?
		bool hasParent(const PackedSeq& seq) const;

		// does this sequence have an extension?
		bool hasChild(const PackedSeq& seq) const;
		
		// print everything
		void printAll() const;
						
		// compute the coordinate of a sequence
		static Coord4 SequenceToCoord4(const Sequence& seq);
		
		// compute the coordinate of a sequence
		static Coord4 SequenceToCoord4(const PackedSeq& pSeq);

	private:
	
		// Compute the transformed (phase space) coordinate of this sequence
		Coord4 SequenceToIndex(const PackedSeq& seq) const;
			
		// Transform a coordinate
		Coord4 CoordToIndex(const Coord4& c) const;
			
		static inline int base2Idx(const char c); 
		
		PhaseSpace();
		
		Bin4D* m_pPhaseSpace;
		Coord4 m_minCoord;
		Coord4 m_maxCoord;
		Coord4 m_size;
		int m_readLength;
		
		bool m_writeEnabled;
};

#endif
