#include "PathDriver.h"

PathDriver::PathDriver(const Sequence& seq, extDirection dir, const PhaseSpace* pPS, const PairRecord* pPR, const SeqRecord* pMR, SeqRecord* pER) : m_seedPath(seq, dir), m_pPhaseSpace(pPS), m_pPairRecord(pPR), m_pMultiplicityRecord(pMR), m_pExtendedSeqRecord(pER)
{
	m_activePaths.push_back(m_seedPath);
	m_pExtendedSeqRecord->addSequence(seq);
}

void PathDriver::addSequence(const Sequence& seq, std::list<Path>& list, bool isSeed)
{

	Path p(seq, m_seedPath.getGrowthDirection());
		
	list.push_back(p);
	m_pExtendedSeqRecord->addSequence(seq);
	//printf("adding %s\n", seq.c_str());
}

void PathDriver::addPairsOfSequence(const Sequence& seq, int position)
{
	return;
	
	// get the pairs of this sequence and add them to the active nodes
	if(m_pPairRecord->checkForPairs(seq))
	{
		const std::vector<Sequence> pairs = m_pPairRecord->getPairs(seq);
		for(const_seq_iter iter = pairs.begin(); iter != pairs.end(); iter++)
		{
			SequencePair p(*iter, position + 200);
			m_validPairs.push_back(p);
		}
	}
}

Path PathDriver::run()
{
	bool stop = false;
	
	while(!stop)
	{
		// extend all the active nodes as far as they can go
		while(extendAllActive(true)) 
		{ 
			trimPairs(m_activePaths.front().getPathLength());
		}
		
		// add all pairs
		printf("all pairs\n");
		for(const_seqPair_iter iter = m_validPairs.begin(); iter != m_validPairs.end(); iter++)
		{
			//printf("pair: %s\n", iter->getSequence().c_str());
		}
	
	
		// merge the paths
		Path seedPath = m_inactivePaths.front();
		return seedPath;
		// store the length before merging to see if it was at all successful
		int premergeLen = seedPath.getPathLength();
		Path subpath(seedPath.getCurrentNode(), seedPath.getGrowthDirection());
		
		std::vector<Path> subpaths = extendSeedPathWithPairs(subpath, 0, 150);
		int bestPathIndex = -1;
		int bestScore = 0;
		int secondBestScore = 0;
		
		for(int i = 0; i < subpaths.size(); i++)
		{
			if(subpaths[i].getPathLength() >= bestScore)
			{
				bestPathIndex = i;
				secondBestScore = bestScore;
				bestScore = subpaths[i].getPathLength();
			}
		}
		
		printf("search returned %d paths best length %d index %d diff %d\n", subpaths.size(), bestScore, bestPathIndex, bestScore - secondBestScore);
		
		if(bestPathIndex == -1 || (bestScore - secondBestScore) < 20)
		{
			return seedPath;
		}
			
		Path bestPath = subpaths[bestPathIndex];
		
		seedPath.mergePath(bestPath, false, false, true);
		
		printf("adding pairs after branch\n");
		if(seedPath.getPathLength() > premergeLen + 30)
		{
			m_activePaths.clear();
			m_inactivePaths.clear();
			//m_validPairs.clear();
			
			int numNodes = seedPath.getPathLength();
			for(int i = numNodes - 200; i < numNodes; i++)
			{
				addPairsOfSequence(seedPath.getNode(i), i);
			}
			
			m_activePaths.push_back(seedPath);
		}
		else
		{
			stop = true;
		}
		printf("done adding pairs after branch\n");
	}

	/*
	for(const_path_iter iter = m_inactivePaths.begin(); iter != m_inactivePaths.end(); iter++)
	{
		printf("full path: %s\n", iter->getSequence().c_str());
	}
	*/
	
	// extension done
	printf("extension step finished %d pairs\n", m_validPairs.size());
	return m_inactivePaths.front();
}

int PathDriver::scoreBranchPath(const Path& path) const
{
	int numNodes = 0;	
}

bool PathDriver::extendAllActive(bool isSeedPath)
{
	// iterate through all the active nodes, extending them one base
	
	std::list<Path>::iterator iter = m_activePaths.begin();
	
	// keep a temporary list of pairs
	std::list<Path> tempList;
	
	while(iter != m_activePaths.end())
	{
		bool isExtended = extendPath(*iter, isSeedPath);

		if(!isExtended)
		{
			//printf("erasing...\n");
			
			// add the path to the inactive nodes
			m_inactivePaths.push_back(*iter);
						
			// could not extend this node, block it
			iter = m_activePaths.erase(iter);
			
			// iter now points to the next element
		}
		else
		{
			// add the pairs of the seed path to the list
			if(isSeedPath)
			{
				addPairsOfSequence(iter->getCurrentNode(), iter->getPathLength());
			}
			
			// increment the iterator
			iter++;
		}	

	}
	
	if(m_activePaths.size() > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

std::vector<Path> PathDriver::extendSeedPathWithPairs(Path seedPath, int distance, int maxDistance)
{
	// conform the valid pairs into a hash for constant time access
	
	// cache this so it doesnt happen every time!
	std::map<Sequence, bool> pairHash;
	for(const_seqPair_iter iter = m_validPairs.begin(); iter != m_validPairs.end(); iter++)
	{
		//printf("adding %s and %s\n", iter->getSequence().c_str(), reverseComplement(iter->getSequence()).c_str());
		pairHash[iter->getSequence()] = true;
		pairHash[reverseComplement(iter->getSequence())] = true;
	}
	
	bool stop = false;
	
	std::vector<Path> retPaths;
	
	while(!stop)
	{
		printf("extending with d %d %s\n", distance, seedPath.getSequence().c_str());
		
		HitRecord hr = m_pPhaseSpace->calculateExtension(seedPath.getCurrentNode(), seedPath.getGrowthDirection());
		if(hr.getNumHits() == 0)
		{
			stop = true;
			retPaths.push_back(seedPath);
		}
		else if(hr.getNumHits() == 1)
		{
			seedPath.addToPath(hr.getFirstHit(), false);
			m_pExtendedSeqRecord->addSequence(hr.getFirstHit());
			distance++;
		}
		else
		{
			
			// check if the ambiguous path is supported by a pair
			int numAmbi = hr.getNumHits();
			printf("num branches %d\n", numAmbi);
			
			
			for(int i = 0; i < numAmbi; i++)
			{
				std::vector<Path> subbranches;
				if(pairHash.find(hr.getHit(i)) != pairHash.end() && !seedPath.contains(hr.getHit(i)))
				{
					Path newPath = seedPath;
					newPath.addToPath(hr.getHit(i), false);
					subbranches = extendSeedPathWithPairs(newPath, distance + 1, maxDistance);
					
					for(int j = 0; j < subbranches.size(); j++)
					{
						retPaths.push_back(subbranches[j]);
					}
				}
				else
				{
					printf("branch not supported %s\n", hr.getHit(i).c_str());	
				}
			}
			
			stop = true;
		}
		
		if(distance > maxDistance)
		{
			stop = true;
			retPaths.push_back(seedPath);
		}	
	}
	
	return retPaths;
}



// the allowInBranch parameter indicates whether we want to allow branching TO a node that has multiple branches incoming
// for seed nodes we don't care about the number of branches going into a node but for pairs
bool PathDriver::extendPath(Path& path, bool allowInBranch)
{
	HitRecord hr = m_pPhaseSpace->calculateExtension(path.getCurrentNode(), path.getGrowthDirection());
	
	//path.print();
	
	if(hr.getNumHits() == 0)
	{
		printf("noext of %s\n", path.getCurrentNode().c_str());
		
		// no ext
		return false;
	}
	else if(hr.getNumHits() == 1)
	{
		path.addToPath(hr.getFirstHit(), false);
		m_pExtendedSeqRecord->addSequence(hr.getFirstHit());
		return true;
	}
	else
	{
		printf("ambiext at %d %s\n", path.getPathLength(), path.getCurrentNode().c_str());
		return false;	
	}
}



void PathDriver::mergePaths()
{
	// try to merge active paths
	std::list<Path> tempPaths;
	
	path_iter iter1 = m_activePaths.begin();
	while(iter1 != m_activePaths.end())
	{
		path_iter iter2 = iter1;
		iter2++;
		
		while(iter2 != m_activePaths.end())
		{
			// check if the sequences are adjacent
			
			bool merged = checkPathsAndMerge(*iter1, *iter2);
			if(merged)
			{
				iter2 = m_activePaths.erase(iter2);
			}
			else
			{
				iter2++;
			}
		}
		iter1++;
	}
	
	// try to merge inactive paths
	iter1 = m_activePaths.begin();
	while(iter1 != m_activePaths.end())
	{
		path_iter iter2 = m_inactivePaths.begin();
		
		while(iter2 != m_inactivePaths.end())
		{
			// check if the sequences are adjacent			
			bool merged = checkPathsAndMerge(*iter1, *iter2);
			if(merged)
			{
				iter2 = m_inactivePaths.erase(iter2);
			}
			else
			{
				iter2++;
			}
		}
		iter1++;
	}	
}

bool PathDriver::checkPathsAndMerge(Path& path1, Path& path2)
{
	SequenceAdjacency rc = checkPathAdjacency(path1, path2);
	if(rc != SA_NONE)
	{
		printf("paths are adjacent:\n");
		path1.print();
		path2.print();
		path1.mergePath(path2, true, false, false);
		return true;
	}
	return false;
}

void PathDriver::trimPairs(int currentPosition)
{
	// remove any pairs that are past the 
	seqPair_iter iter = m_validPairs.begin();
	
	while(iter != m_validPairs.end()) 
	{
		if(iter->getMaxPos() < currentPosition)
		{
			iter = m_validPairs.erase(iter);
		}
		else
		{
			iter++;
		}
	}
}

SequenceAdjacency PathDriver::checkPathAdjacency(const Path& path1, const Path& path2) const
{	
	
	SequenceAdjacency rc = SA_NONE;
	
	// check the sense direction of path 1 against antisense of path2
	const Sequence& sSeq1 = path1.getCurrentNode();
	const Sequence& sSeq2 = path2.getCurrentNode();	
	
	if(checkSequenceAdjacency(sSeq1, sSeq1, SENSE))
	{
		rc = SA_SAME_SENSE;
	}
	
	if(checkSequenceAdjacency(sSeq1, reverseComplement(sSeq2), SENSE))
	{
		rc = SA_RC_SENSE;
	}
		
	return rc;
}

bool PathDriver::checkSequenceAdjacency(const Sequence& seq1, const Sequence& seq2, extDirection dir) const
{
	// chop the sequence by 1 base in the direction of the extension
	int commonLen = seq1.length() - 1;
	
	Sequence commonSeq1;
	Sequence commonSeq2;
	
	if(dir == SENSE)
	{
		commonSeq1 = seq1.substr(1,commonLen);
		commonSeq2 = seq2.substr(0, commonLen);
	}
	else
	{
		commonSeq1 = seq1.substr(0,commonLen);
		commonSeq2 = seq2.substr(1, commonLen);
	}
	
	if(commonSeq1 == commonSeq2)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void PathDriver::printAll(Writer& writer) const
{
	for(const_path_iter iter = m_activePaths.begin(); iter != m_activePaths.end(); iter++)
	{
		writer.writeContig(iter->getSequence().c_str());	
	}

	for(const_path_iter iter = m_inactivePaths.begin(); iter != m_inactivePaths.end(); iter++)
	{
		writer.writeContig(iter->getSequence().c_str());	
	}
}
