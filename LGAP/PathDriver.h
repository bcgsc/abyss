#ifndef PATHDRIVER_H
#define PATHDRIVER_H

#include <list>
#include "Path.h"
#include "HitRecord.h"
#include "PhaseSpace.h"
#include "PairRecord.h"
#include "Writer.h"
#include "SeqRecord.h"
#include "SequencePair.h"

typedef std::list<Path>::iterator path_iter;
typedef std::list<Path>::const_iterator const_path_iter;

typedef std::list<SequencePair>::iterator seqPair_iter;
typedef std::list<SequencePair>::const_iterator const_seqPair_iter;

class PathDriver
{
	public:
	
		// constructor
		PathDriver(const Sequence& seq, extDirection dir, const PhaseSpace* pPS, const PairRecord* pPR, const SeqRecord* pMR, SeqRecord* pER);
		
		// add a sequence to extend
		void addSequence(const Sequence& seq, std::list<Path>& list, bool isSeed);
		
		// add the pairs of the sequence to extend
		void addPairsOfSequence(const Sequence& seq, int position);
		
		// extend all the nodes as far as they can go (until they hit an ambiguous or already seen sequence)
		bool extendAllActive(bool isSeedPath);
		
		// extend a path in a single direction
		bool extendPath(Path& path, bool allowInBranch);
		
		// extend the seed path using paired info
		std::vector<Path> extendSeedPathWithPairs(Path seedPath, int distance, int maxDistance);
		
		// extend paths given pair support
		std::vector<Path> extendPathSupported(Path currPath, extDirection dir);
		
		// check the adjacency between two sequences
		SequenceAdjacency checkPathAdjacency(const Path& path1, const Path& seq2) const;
		
		bool checkSequenceAdjacency(const Sequence& seq1, const Sequence& seq2, extDirection dir) const;
		
		// check and merge
		bool checkPathsAndMerge(Path& path1, Path& path2);
		
		
		// run one step of the extension
		Path run();
		
		// print all the current paths
		void printAll(Writer& writer) const;
		
		// merge paths
		void mergePaths();
		
		// trim paths
		void trimPairs(int currentPosition);
		
		// score branch path
		int scoreBranchPath(const Path& path) const;
		
	private:
	
		// the initial path to grow from
		Path m_seedPath;
		
		// the main record of all the sections that can be seen
		const PhaseSpace* m_pPhaseSpace;
		
		// the main sequence->pair mapping structure
		const PairRecord* m_pPairRecord;
		
		// the multiplicity of each sequence
		const SeqRecord* m_pMultiplicityRecord;
		
		// the record containing which sequences have already been extended
		SeqRecord* m_pExtendedSeqRecord;
		
		// active/inactive paths
		std::list<Path> m_activePaths;
		std::list<Path> m_inactivePaths;

		// a list of the currently valid pairs
		std::list<SequencePair> m_validPairs;

};

#endif
