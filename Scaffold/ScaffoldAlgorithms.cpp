#include "ScaffoldAlgorithms.h"

namespace ScaffoldAlgorithms
{
	
//
// generate the master start list
//
void generateStartList(ISequenceCollection* seqCollection, std::vector<ContigStart>& startList)
{
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		// a contig (single-end case) can start by either:
		// 1) a sequence endpoint where one half of the sequence has no extension
		// 2) an ambiguious extension on the opposite direction
				
		// check case 1
		extDirection dir;
		SeqContiguity status = AssemblyAlgorithms::checkSeqContiguity(seqCollection, *iter, dir);
		if(status != SC_INVALID)
		{
			if(status == SC_ISLAND)
			{
				ContigStart newStart;
				newStart.seq = *iter;
				newStart.state = CSS_ISLAND;
				newStart.dir = SENSE; // meaningless since it is an island so set default
				startList.push_back(newStart);				
				continue;
			}
			else if(status == SC_ENDPOINT)
			{
				// the sequence is discontinus
				ContigStart newStart;
				newStart.seq = *iter;
				newStart.state = CSS_ENDPOINT;
				newStart.dir = dir;
				startList.push_back(newStart);
				
				continue;
			}
			
			if(status == SC_CONTIGUOUS)
			{
				ExtensionRecord extRec;
				int multiplicity;
				
				// not network friendly
				bool success = seqCollection->getSeqData(*iter, extRec, multiplicity);
				assert(success);
				(void)success;
				
				
				for(int i = 0; i <= 1; ++i)
				{
					extDirection dir = (i == 0) ? SENSE : ANTISENSE;
					
					if(extRec.dir[dir].IsAmbiguous())
					{
						// Add ambiguous starts
						ContigStart newStart;
						newStart.seq = *iter;
						newStart.state = CSS_AMBIGUOUS;
						newStart.dir = !dir; // opposite direction
						//startList.push_back(newStart);				
						
						// add a start for each ambiguous seq
						PSequenceVector newSeqs;
						AssemblyAlgorithms::generateSequencesFromExtension(*iter, dir, extRec.dir[dir], newSeqs);
						
						for(PSequenceVector::iterator nsIter = newSeqs.begin(); nsIter != newSeqs.end(); ++nsIter)
						{
							ContigStart newStart;
							newStart.seq = *nsIter;
							newStart.state = CSS_AMBIGUOUS;
							newStart.dir = dir;
							startList.push_back(newStart);
						}									
					}
				}
			}			
		}
	}
}

//
//
//
bool processNonlinearExtensionForBranch(ISequenceCollection* seqCollection, PairRecord* pPairRecord, BranchRecord& branch, PackedSeq& currSeq, ExtensionRecord extensions, int multiplicity)
{
	extDirection dir = branch.getDirection();
	//extDirection oppDir = oppositeDirection(dir);
	
	if(branch.doLengthCheck() && ((int)branch.getLength() > (int)branch.getMaxLength())) // Check if the branch has extended past the max trim length
	{
		// Stop the branch
		branch.terminate(BS_TOO_LONG);
	}
	else if(branch.hasLoop())
	{
		branch.terminate(BS_LOOP);
	}
	/*
	else if(extensions.dir[oppDir].IsAmbiguous()) // Does this sequence split TO the former node?
	{
		//printf("stopped because of reverse branch\n");
		// There is a reverse ambiguity to this branch, stop the branch without adding the current sequence to it
		branch.terminate(BS_AMBI_OPP);
	}*/
	else if(!extensions.dir[dir].HasExtension()) 
	{
		// no extenstion, add the current sequence and terminate the branch
		branch.addSequence(currSeq);
		branch.setMultiplicity(currSeq, multiplicity);
		branch.terminate(BS_NOEXT);
	}
	else if(extensions.dir[dir].IsAmbiguous())
	{
		branch.addSequence(currSeq);
		branch.setMultiplicity(currSeq, multiplicity);
					
		// try to resolve this spot of ambiguity
		PackedSeq chosenNode;
		bool resolved = deconvolvePaths(seqCollection, pPairRecord, branch, currSeq, dir, chosenNode);
		if(resolved)
		{
			currSeq = chosenNode;
		}
		else
		{
			// terminate
			branch.terminate(BS_AMBI_SAME);
		}
	}
	else
	{
		// Add the sequence to the branch
		branch.addSequence(currSeq);
		branch.setMultiplicity(currSeq, multiplicity);
		
		// generate the new current sequence from the extension
		//printf("currseq: %s ", currSeq.decode().c_str());
		PSequenceVector newSeqs;

		AssemblyAlgorithms::generateSequencesFromExtension(currSeq, dir, extensions.dir[dir], newSeqs);
		assert(newSeqs.size() == 1);
		currSeq = newSeqs.front();
		//printf("newseq: %s \n", currSeq.decode().c_str());
	}
	
	return branch.isActive();
}

bool deconvolvePaths(ISequenceCollection* seqCollection, PairRecord* pPairRecord, const BranchRecord& currBranch, const PackedSeq& branchPoint, extDirection dir, PackedSeq& chosenNode)
{
	// Build the tree from this endmer
	const int insertLen = 250;
	int lookForwardDepth = insertLen / 2;
	int lookBackDepth = insertLen;

	ParentTree parentTree(branchPoint, dir, seqCollection, lookForwardDepth);

	// Get pairs
	int numNodesInBranch = currBranch.getLength();
	size_t firstIndex;
	size_t lastIndex = numNodesInBranch - 1;
							
	if((int)lastIndex - lookBackDepth < 0)
	{
		firstIndex = 0;
	}
	else
	{
		firstIndex = lastIndex - lookBackDepth;
	}
									
	PSeqSet children = parentTree.getRootChildren();					
	const double minSupport = 10.0f;
	
	PSequenceVector allPairs;
	
	int numGood = 0;
	
	for(size_t idx = firstIndex; idx < lastIndex; ++idx)
	{		
		double weight = 1.0f;
		PSequenceVector seqPairs = pPairRecord->getPairs(currBranch.getSeqByIndex(idx));
		
		for(PSequenceVector::iterator pairIter = seqPairs.begin(); pairIter != seqPairs.end(); pairIter++)
		{
			// Check if this pair is supporting multiple branches
			// If it is, do not add the pairs to the set
			parentTree.scoreRootChildrenOnly(*pairIter, weight);
			//parentTree.addScore(*pairIter, weight);
		}
		
		
		// Check the number of branches supported by the pairs of this seq
		int numSupported = 0;
		int suppIdx = 0;
		int currIdx = 0;
		for(PSeqSet::iterator testIter = children.begin(); testIter != children.end(); ++testIter)
		{
			double currScore = parentTree.getScore(*testIter).subtreeScore;	
			if(currScore > 0)
			{
				suppIdx = currIdx;
				numSupported++;
			}
			currIdx++;
		}					

		if(numSupported == 1)
		{
			for(PSequenceVector::iterator pairIter = seqPairs.begin(); pairIter != seqPairs.end(); pairIter++)
			{
				allPairs.push_back(*pairIter);
			}
			numGood++;
		}
		parentTree.clearRootChildrenScores();
		//parentTree.resetScore();
	}
	
	printf("%d of %d sequences are usable\n", numGood, lastIndex - firstIndex);
	
	int count = 0;
	for(PSequenceVector::iterator pairIter = allPairs.begin(); pairIter != allPairs.end(); pairIter++)
	{
		count++;
		parentTree.scoreRootChildrenOnly(*pairIter, 1.0f);
		//parentTree.addScore(*pairIter, 1.0f);
	}
	
	parentTree.print(3);
	
	parentTree.printRootChildren();

	PSeqSet::iterator bestIter;
	
	int numSupported = 0;
	for(PSeqSet::iterator testIter = children.begin(); testIter != children.end(); ++testIter)
	{
		double currScore = parentTree.getScore(*testIter).subtreeScore;	
		if(currScore > minSupport)
		{
			bestIter = testIter;
			numSupported++;
		}
	}
	
	if(numSupported == 1)
	{
		printf("path resolved\n");
		chosenNode = *bestIter;
		return true;
	}
	else
	{
		printf("path not resolved\n");
		return false;	
	}
}

};
