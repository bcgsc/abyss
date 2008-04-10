#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include <iostream>
#include <fstream>
#include <map>
#include "CommonDefs.h"
#include "ISequenceCollection.h"
#include "AssemblyAlgorithms.h"

// CORIEN_SAME = the contigs are from the same strand
// CORIEN_OPP = the contigs are from different strands
// CORIEN_AMBI = the orientation is ambiguous
enum ContigOrientation
{
	CORIEN_SAME,
	CORIEN_OPP,
	CORIEN_AMBI
};

struct range
{
	int start;
	int end;	
};

// ContigPlacement is with respect to the first contig (reference)
// C1---C2 -> CORDER_LEFT
// C2---C1 -> CORDER_RIGHT
 
enum ContigOrder
{
	CORDER_LEFT,
	CORDER_RIGHT	
};

enum LinkType
{
	LT_STRONG,
	LT_WEAK,
	LT_REMOVED
};

struct ReadAlign
{
	ReadID id;
	ContigID contig;
	int pos;
	bool isRC;
};

struct PairAlign
{
	ReadAlign pairs[2];
	bool invalid;
};

struct ContigLinkage
{
	ContigID masterID;
	ContigID slaveID;
	ContigOrientation orientation;
	ContigOrder order;
	int distance;
	int numPairs;
	LinkType type;
	
	// Flag indicating there is no good linkage
	bool noLink;
};
	

// Typedefs
typedef std::vector<PairAlign> PairAlignVec;
typedef std::vector<ReadAlign> AlignVec;

typedef std::map<ReadID, AlignVec> AlignmentMap;
typedef std::map<ReadID, ReadID> PairingMap;
typedef std::set<ReadID> ReadSet;
typedef std::map<ContigID, ReadSet> ContigReadMap;
typedef std::vector<ContigLinkage> LinkVec;
typedef std::vector<Sequence> SeqVec;

typedef std::map<ContigID, PairAlignVec> ContigPairVecMap;

typedef std::vector<double> PDF;
typedef std::map<int, int> histogram;



// Iterators
typedef AlignmentMap::iterator AMIter;
typedef PairingMap::iterator PMIter;
typedef ReadSet::iterator RSIter;
typedef ContigReadMap::iterator CRMIter;
typedef ContigPairVecMap::iterator CPVMIter;
typedef PairAlignVec::iterator PAVIter;
typedef LinkVec::iterator LinkIter;
typedef SeqVec::iterator SeqVecIter;
typedef AlignVec::iterator AVIter;


typedef AlignmentMap::const_iterator ConstAMIter;
typedef PairingMap::const_iterator ConstPMIter;
typedef ReadSet::const_iterator ConstRSIter;
typedef ContigReadMap::const_iterator ConstCRMIter;

class Scaffold
{
	public:
		Scaffold(std::string readsFile, std::string contigFile, int readLen, int kmer, bool bReadAlignments, std::string alignmentsFile);
		
		//IO Functions
		void ReadSequenceReads(std::string file);
		void ReadPairs(std::string file);
		void ReadAlignments(std::string file);
		void ReadContigs(std::string file);
		
		// Generate the empirical distribution of pair distances
		void GenerateEmpDistribution();
		void ConvertHistToPDF(const histogram& h);
		
		// Determine the adjacency information between contigs
		bool AttemptMerge(ContigID contigID);
		ContigLinkage GenerateLinkage(ContigID contigID0, ContigID contigID1, PairAlignVec& paVec);

		
	private:
	
		// Generate alignments of the read set against the contigs
		void GenerateAlignments(PSequenceVector& seqs, ContigMap& contigs);
		 
		// Refine the linkages, attempting to upgrade weak links to strong links
		bool RefineLinkages(LinkVec& links);
	
		// Generate a graphwiz graph of the linkages around the particular contig
		void GenerateGraph(ContigID contigID);
		
		// Merge contigs
		int Merge(Sequence& leftContig, Sequence& rightContig, int distance, Sequence& merged);
		
		// Sub assemble the paired reads
		SeqVec SubAssemble(PSequenceVector& seqs, Sequence startNode, int maxDistance);
		
		// Recursively assemble
		SeqVec AssembleRecursive(ISequenceCollection* pSC, extDirection dir, PackedSeq start, int d, int maxDistance);
		
		// Check if the links are consistent with the chosen best link
		bool CheckConsistency(ContigLinkage bestLink, LinkVec& alllinks);

		// Align contigs using the input position as a guess to the alignment
		int alignContigs(const Sequence& leftContig, const Sequence& rightContig, int guess, int range, int& retScore);
		
		// Get all the reads of the pairs that are of the specified complement
		void GetEndPairs(ContigID contigID, bool compPairs, PSequenceVector& outSeqs);
		
		// Update the positions of the reads on the master contig
		void UpdateMasterReads(ContigID contigID, int offset, const Sequence& origSeq, const Sequence& merged);

		// Update the positions of the reads on the slave contig
		void UpdateSlaveReads(ContigID slaveID, ContigID masterID, int offset, bool isFlipped, const Sequence& origSeq, const Sequence& merged);
		
		// Realign the pairs of the reads on the specified contig
		void RealignContigPairs(ContigID contigID);
	
		// Determine the orientation between contigs
		ContigOrientation DetermineOrientation(PairAlignVec& contigPairs);

		// Determine the order of the contigs
		ContigOrder DetermineOrder(PairAlignVec& contigPairs);
		
		// Reverse the positions of the pairs for the second contig pairs
		void ReverseSecondContigPairs(PairAlignVec& contigPairs, int contigLength);
		
		// Estimate the distance between the contigs
		int EstimateDistanceBetweenContigs(PairAlignVec& contigPairs, ContigOrder order, Sequence& contig1, Sequence& contig2);
		
		// Get all the pairs between the specific contig and any other contigs and places them in cpvMap
		void GenerateUniquePairAlignments(ContigID contigID, ContigPairVecMap& cpvMap);
		
		// Populate the pairAlign data structure if the pairs are unique and both are aligned
		// Returns true if unique/aligned, false otherwise
		bool GetUniquePairAlign(ReadID readID, PairAlign& pairAlign);
		
		// Maximum Likelihood Estimator functions
		int MaxLikelihoodEst(std::vector<int>& pairDistance, PDF& pdf);
		
		// Compute the likelihood of the distribution
		double ComputeLikelihood(int d, std::vector<int>& testDist, PDF& pdf);
		
		// Get the alignments for a particular read
		AlignVec GetAlignmentsForRead(ReadID id);
		
		// Get the sequence for a particular read
		PackedSeq GetSequenceForRead(ReadID id);		
		
		// Get the ID of the pair of the read
		ReadID GetPairID(ReadID id);
		
		// Build an alignment structure
		ReadAlign BuildReadAlign(ReadID id, std::string contig, int position, bool isRC);
		
		// Print a read alignment
		void PrintReadAlign(ReadAlign& ra);
		
		// Print a linkage
		void PrintLinkage(ContigLinkage& link);
		
		// Calculate the amount of overlap between the ranges
		int OverlapRanges(const range r1, const range r2);
		
		// Generate the standard deviation of the estimate
		double GetStdDevOfEstimate(int n);
		
		// generate the maximum coordinate set for the specified number if deviations
		range GenerateRange(int distance, int size, int n, int numDevs);
		
		// Print a gviz node to the file handle
		void OutputGVizNode(std::ofstream& ostr, ContigLinkage& link);
		
		// Write all the alignments to a file
		void WriteAlignments(std::string filename);
		
		// All the mapping datastructures needed
		AlignmentMap m_alignMap;
		//PairingMap m_pairMap; 
		ContigReadMap m_contigReadMap;
		ContigMap m_contigMap;
		PSequenceVector m_readVec;
		
		PDF m_pdf;
		double m_stdDev;
		
		int m_readLen;
		int m_kmer;
		
		static const int STRONG_LINK_CUTOFF = 10;
		static const int SUB_ASSEMBLY_K = 12;
};


int CompareLinkagesByDistance(const ContigLinkage& l1, const ContigLinkage& l2);
int CompareLinkagesByDistanceDesc(const ContigLinkage& l1, const ContigLinkage& l2);

#endif
