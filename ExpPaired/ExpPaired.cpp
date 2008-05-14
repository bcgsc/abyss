#include <stdio.h>

#include <vector>
#include <stdio.h>
#include <deque>
#include <iostream>
#include <fstream>
#include "ExpPaired.h"
#include "CommonUtils.h"
#include "DotWriter.h"
#include "ISequenceCollection.h"
#include "Options.h"
#include "SequenceCollection.h"
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"
#include "Timer.h"
#include "PairRecord.h"
#include "NodeScores.h"
#include "ParentTree.h"
#include "PackedSeqReader.h"

int main(int /*argc*/, char* const* argv)
{	
	int kmerSize = atoi(argv[1]);
	std::string graphPacked = argv[2];
	std::string pairsFile = argv[3];
	
	printf("kmer %d graph: %s pairs: %s\n", kmerSize, graphPacked.c_str(), pairsFile.c_str());
	
	// Create the sequence hash
	SequenceCollectionHash* pSC = new SequenceCollectionHash();
	
	loadSeqHash(graphPacked, pSC);
	
	printf("finalizing\n");
	pSC->finalize();

	AssemblyAlgorithms::generateAdjacency(pSC);
	
	printf("loading pairs\n");
	PairRecord pairRec;
	loadPairs(pairRec, pairsFile, kmerSize);
	
	FastaWriter* originalWriter = new FastaWriter("original_contigs.fa");
	FastaWriter* newWriter = new FastaWriter("new_contigs.fa");
		
	printf("assembling\n");
	//while(assemble_paired(pSC, pairRec, originalWriter, newWriter) > 0);
	
	assemble2(pSC, pairRec, originalWriter, newWriter);
	
	delete originalWriter;
	delete newWriter;	
}

void loadSeqHash(std::string packedSeqFile, ISequenceCollection* pSS)
{
	printf("Reading packed seq binary file %s\n", packedSeqFile.c_str());
	PackedSeqReader* reader = new PackedSeqReader(packedSeqFile.c_str());

	bool stop = false;
	while(!stop)
	{
		PSequenceVector seqs;
		stop = !reader->ReadSequences(seqs);
		for(PSequenceVectorIterator iter = seqs.begin(); iter != seqs.end(); iter++)
		{			
			pSS->add(*iter);
		}
	}
		
		delete reader;
}

void loadPairs(PairRecord& pairRec, std::string pairFile, int kmerSize)
{
	Timer timer("Read Pairs");
	FastaReader* reader = new FastaReader(pairFile.c_str());
	int count = 0;
	while(reader->isGood())
	{
		Sequence seq1 = reader->ReadSequence();
		Sequence seq2 = reader->ReadSequence();
		
		// break into kmers
		assert(seq1.length() == seq2.length());
		
		int len = seq1.length();
		assert(kmerSize <= len);
		
		for(int i = 0; i < len - kmerSize  + 1; i++)
		{
			PackedSeq sub1(seq1.substr(i, kmerSize));
			PackedSeq sub2(seq2.substr(i, kmerSize));
			
			pairRec.addPairs(sub1, sub2);
			count++;
		}
	}
	printf("Read %d pairs\n", count);
	delete reader;
}

int assemble_paired(ISequenceCollection* seqCollection, PairRecord& pairRecord, FastaWriter* originalWriter, FastaWriter* newWriter)
{

	int numAssembled = 0;
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	int id = 0;
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		bool first = true;
		extDirection dir;
		// dir will be set to the trimming direction if the sequence can be trimmed
		SeqContiguity status = AssemblyAlgorithms::checkSeqContiguity(seqCollection, *iter, dir);

		if(status == SC_INVALID || status == SC_CONTIGUOUS)
		{
			continue;
		}
		else if(status == SC_ISLAND)
		{
			// singleton, ignore for now
			continue;
		}
		
		// The sequence is an endpoint, begin extending it
		// Passing -1 into the branch will disable the length check
		BranchRecord currBranch(dir, -1);
					
		PackedSeq currSeq = *iter;
		
		while(currBranch.isActive())
		{		
			// Get the extensions for this sequence, this function populates the extRecord structure
			ExtensionRecord extRec;
			bool success = seqCollection->getExtensions(currSeq, extRec);
			assert(success);
			
			// process the extension record and extend the current branch, this function updates currSeq on successful extension
			bool inBranchDetected = false;
			processExtensionPaired(currBranch, currSeq, extRec, inBranchDetected);
			
			/*
			// Branch has hit ambiguitity
			if(currBranch.getState() == BS_AMBI_OPP)
			{
				currBranch.setState(BS_ACTIVE);
			}
			*/
			if(currBranch.getState() == BS_AMBI_SAME)
			{
				printf("Ambiguous branch found at distance %zu, trying to resolve\n", currBranch.getLength());
				printf("Current seq %s\n", currSeq.decode().c_str());
				Sequence seq;
				currBranch.buildContig(seq);
				
				if(first)
				{
					originalWriter->WriteSequence(seq, id, 0.0f);
					first = false;
				}
				
				printf("Current contig: %s\n", seq.c_str());
				//const int maxDepth = 200;
				const int lookBackDepth = 200;
				const int lookForwardDepth = 10;
				ParentTree parentTree(currSeq, dir, seqCollection, lookForwardDepth);
				
				// Get all the pairs of the sequences leading up to the current sequence
								
				PairScoreVec allPairs;
				
				int numNodesInBranch = currBranch.getLength();
				int firstIndex = numNodesInBranch - lookBackDepth;
				int lastIndex = numNodesInBranch - 100;
				if(firstIndex < 0)
				{
					firstIndex = 0;
				}
				
				for(int idx = firstIndex; idx < lastIndex; ++idx)
				{
					PackedSeq firstHalf = currBranch.getSeqByIndex(idx);
					
					/*
					int mult = seqCollection->getMultiplicity(firstHalf);
					mult = mult * mult;
					printf("firsthalf %s, mult %d\n", firstHalf.decode().c_str(), mult);
					*/
					double weight = 1.0f;
					PSequenceVector seqPairs = pairRecord.getPairs(firstHalf);
					for(PSequenceVector::iterator pairIter = seqPairs.begin(); pairIter != seqPairs.end(); pairIter++)
					{
						PairScore ps;
						ps.weight = weight;
						ps.seq = *pairIter;
						
						allPairs.push_back(ps);
					}
				}
				

				// score all the pairs
				printf("scoring %zu pairs\n", allPairs.size());
				

				
				int count = 0;
				for(PairScoreVec::iterator pairIter = allPairs.begin(); pairIter != allPairs.end(); pairIter++)
				{
					count++;
					parentTree.addScore(pairIter->seq, pairIter->weight);
				}
				
				
				parentTree.print(10);
				
				PSeqSet children = parentTree.getRootChildren();
				parentTree.printRootChildren();
				PSeqSet::iterator bestIter;
								
				const double minSupport = 10.0f;
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
				// choose a branch
				
				//const double cutoff = 0.55f;
				//double frac = bestScore / sumScore;
				//if(children.size() > 0 && frac > cutoff)
				if(numSupported == 1)
				{
					printf("Reactivating branch\n");
					//currBranch.addSequence(children[bestIndex]);
					currBranch.setState(BS_ACTIVE);
					currSeq = *bestIter;
				}
				else
				{
					printf("Branch finished for real\n");
					currBranch.setState(BS_AMBI_SAME);
					Sequence seq;
					currBranch.buildContig(seq);
					printf("last contig: %s\n", seq.c_str());					
				}
			}
		}

				
				
		Sequence seq;
		
		// If the assembly stopped because of ambiguity, split the last sequence
		if(currBranch.getState() == BS_AMBI_SAME)
		{
			AssemblyAlgorithms::removeExtensionsToSequence(seqCollection, currBranch.getLastSeq(), currBranch.getDirection());
			seqCollection->clearExtensions(currBranch.getLastSeq(), currBranch.getDirection());	
		}
		
		//hack
		currBranch.terminate(BS_NOEXT);
		AssemblyAlgorithms::processTerminatedBranchAssemble(seqCollection, currBranch, seq);
		
		newWriter->WriteSequence(seq, id, 0.0f);
		
		numAssembled++;
		id++;
	}
	
	return numAssembled;
}

int assemble2(ISequenceCollection* seqCollection, PairRecord& /*pairRecord*/, FastaWriter* /*originalWriter*/, FastaWriter* newWriter)
{

	int numAssembled = 0;

	// Build the list of contig starts
	startList contigStartList;
	
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{	
		extDirection dir;
		SeqContiguity status = AssemblyAlgorithms::checkSeqContiguity(seqCollection, *iter, dir);
		if(status == SC_ENDPOINT)
		{
			ContigStart start;
			start.dir = dir;
			start.seq = *iter;
			contigStartList.push_back(start);			
		}
	}
		
	int id = 0;
	
	while(!contigStartList.empty())
	{
		printf("Num starts: %zu\n", contigStartList.size());
		startList newStarts;
		
		for(startList::iterator iter = contigStartList.begin(); iter != contigStartList.end(); ++iter)
		{
			BranchRecord currBranch(iter->dir, -1);
			PackedSeq currSeq = iter->seq;
			extDirection dir = iter->dir;

			SeqContiguity status = AssemblyAlgorithms::checkSeqContiguity(seqCollection, currSeq, dir);
			if(status == SC_INVALID)
			{
				continue;
			}
				
			while(currBranch.isActive())
			{		
				// Get the extensions for this sequence, this function populates the extRecord structure
				ExtensionRecord extRec;
				bool success = seqCollection->getExtensions(currSeq, extRec);
				assert(success);
				
				// process the extension record and extend the current branch, this function updates currSeq on successful extension
				bool inBranchDetected = false;
				processExtensionPaired(currBranch, currSeq, extRec, inBranchDetected);
				
				if(currBranch.getState() == BS_AMBI_SAME)
				{
					const int lookForwardDepth = 10;
					ParentTree parentTree(currSeq, dir, seqCollection, lookForwardDepth);

					
					PSeqSet children = parentTree.getRootChildren();
					parentTree.printRootChildren();
					PSeqSet::iterator bestIter;

					for(PSeqSet::iterator testIter = children.begin(); testIter != children.end(); ++testIter)
					{
						ContigStart start;
						start.dir = currBranch.getDirection();
						start.seq = *testIter;
						newStarts.push_back(start);
					}
				}
			}
			
			

			//hack
			currBranch.terminate(BS_NOEXT);
			Sequence seq;
			AssemblyAlgorithms::processTerminatedBranchAssemble(seqCollection, currBranch, seq);
			
			newWriter->WriteSequence(seq, id, 0.0f);
			
			numAssembled++;
			id++;						
		}
		
		contigStartList = newStarts;	
	}
	return numAssembled;
}

//
BranchState processExtensionPaired(BranchRecord& branch, PackedSeq& currSeq, ExtensionRecord extensions, bool inBranchDetected)
{
	extDirection dir = branch.getDirection();
	extDirection oppDir = oppositeDirection(dir);
	
	
	// Check if there is an inbranch
	if(extensions.dir[oppDir].IsAmbiguous())
	{
		inBranchDetected = true;	
	}
	
	if(branch.hasLoop())
	{
		branch.terminate(BS_LOOP);
		return BS_LOOP;
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
		branch.terminate(BS_NOEXT);
		return BS_NOEXT;
	}
	else if(extensions.dir[dir].IsAmbiguous())
	{
		// this branch has an ambiguous extension, add the current sequence and terminate
		branch.addSequence(currSeq);
		branch.terminate(BS_AMBI_SAME);
		return BS_AMBI_SAME;
	}
	else
	{
		// Add the sequence to the branch
		branch.addSequence(currSeq);
		
		// generate the new current sequence from the extension
		//printf("currseq: %s ", currSeq.decode().c_str());
		PSequenceVector newSeqs;

		AssemblyAlgorithms::generateSequencesFromExtension(currSeq, dir, extensions.dir[dir], newSeqs);
		assert(newSeqs.size() == 1);
		currSeq = newSeqs.front();
	}
	
	return BS_ACTIVE;
}
