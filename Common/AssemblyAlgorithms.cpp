#include <iostream>
#include "AssemblyAlgorithms.h"
#include "Options.h"
#include "SequenceCollectionHash.h"
#include "Timer.h"

namespace AssemblyAlgorithms
{

//
// Generate a hitrecord for a sequence
//
HitRecord calculateExtension(ISequenceCollection* seqCollection,
		const PackedSeq& currSeq, extDirection dir)
{
	
	// Create the return structure
	HitRecord hitRecord;
	
	// Check for the existance of the 4 possible extensions
	for(int i  = 0; i < NUM_BASES; i++)
	{
		char currBase = BASES[i];
		
		ResultPair hasExt = seqCollection->checkExtension(currSeq, dir, currBase);
		//if(printAll) printf("%s has extension(%d) to: %c ? %d\n", currSeq.decode().c_str(), dir, currBase, hasExt.forward || hasExt.reverse);
		
		// Does this sequence have an extension?
		if(hasExt.forward || hasExt.reverse)
		{
			PackedSeq extSeq(currSeq);
			extSeq.rotate(dir, currBase);
			
			// is there a forward extension?
			if(hasExt.forward)
			{
				hitRecord.addHit(extSeq, false);
			}
			else
			{
				// extension is of the reverse complement
				hitRecord.addHit(extSeq, true);	
			}
		}
	}
		
	return hitRecord;
}

//
// Generate a hitrecord for a sequence
//
void generateSequencesFromExtension(const PackedSeq& currSeq, extDirection dir, SeqExt extension, PSequenceVector& outseqs)
{
	
	// Create the return structure
	PSequenceVector extensions;

	PackedSeq extSeq(currSeq);
	extSeq.rotate(dir, 'A');
	
	// Check for the existance of the 4 possible extensions
	for(int i  = 0; i < NUM_BASES; i++)
	{
		char currBase = BASES[i];
	
		// Does this sequence have an extension?
		if(extension.CheckBase(currBase))
		{
			extSeq.setLastBase(dir, currBase);
			outseqs.push_back(extSeq);
		}
	}
}

//
// Function to load sequences into the collection
//
void loadSequences(ISequenceCollection* seqCollection,
		std::string fastaFile, int /*readLength*/, int kmerSize)
{
	Timer timer("LoadSequences");
	FastaReader* reader = new FastaReader(fastaFile.c_str());
	
	int count = 0;
	
	// Read the sequences and add them to the network sequence space
	int64_t idNum = 0;
	size_t lastNum = 0;
	while(reader->isGood())
	{
		Sequence seq = reader->ReadSequence();
		int len = seq.length();
		assert(kmerSize <= len);
		//assert(seq.getSequenceLength() == pairSeq.getSequenceLength());
		
		//int l = pairSeq.getSequenceLength();

		for(int i = 0; i < len - kmerSize  + 1; i++)
		{
			PackedSeq sub(seq.substr(i, kmerSize));
			
			//PackedSeq pairSub = pairSeq.subseq(l - kmerSize - i, kmerSize);

			// Add the sequence to the sequence collection
			seqCollection->add(sub);

		}
	
		if(count % 100000 == 0)
		{
			size_t numseqs = seqCollection->count();
			std::cout << "read " << count << " sequences (" << seqCollection->count() << " delta: " << numseqs - lastNum << ")" << std::endl;
			std::cout.flush();
			lastNum = numseqs;
		}
		count++;
			
		idNum++;
		seqCollection->pumpNetwork();
	}
	
	delete reader;
	reader = 0;
}

//
// Generate the adjacency information for each sequence in the collection
//
void generateAdjacency(ISequenceCollection* seqCollection)
{
	Timer timer("GenerateAdjacency");
	printf("generating adjacency info\n");
	

	int count = 0;
	int numBasesSet = 0;
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		if(count % 1000000 == 0)
		{
			printf("generated for %d\n", count);
		}
		count++;
		
		//printf("gen for: %s\n", iter->decode().c_str());
		for(int i = 0; i <= 1; i++)
		{
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			SeqExt extension;
			PackedSeq testSeq(*iter);
			testSeq.rotate(dir, 'A');
			for(int j = 0; j < NUM_BASES; j++)
			{
				char currBase = BASES[j];
				testSeq.setLastBase(dir, currBase);
				
				if(seqCollection->exists(testSeq))
				{		
					extension.SetBase(currBase);
					numBasesSet++;
				}
			}
			seqCollection->setExtension(*iter, dir, extension);	
		}
		
		//iter->printExtension();
		seqCollection->pumpNetwork();
	}
	
	printf("adjacency set %d bases\n", numBasesSet);
}

//
//
//
void splitAmbiguous(ISequenceCollection* seqCollection)
{
	Timer timer("SplitAmbiguous");
	printf("splitting ambiguous reads\n");
	int count = 0;
	int numSplit = 0;
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		if(iter->isFlagSet(SF_DELETE))
		{
			continue;
		}
		
		if(count % 1000000 == 0)
		{
			printf("split %d\n", count);
		}
		count++;
		
		for(int i = 0; i <= 1; i++)
		{
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			HitRecord hr = calculateExtension(seqCollection, *iter, dir);
			if(hr.getNumHits() > 1)
			{
				removeExtensionsToSequence(seqCollection, *iter, dir);
				seqCollection->clearExtensions(*iter, dir);	
				numSplit++;			
			}
		}
		seqCollection->pumpNetwork();
	}
	printf("Split %d ambiguous nodes\n", numSplit);
}

int popBubbles(ISequenceCollection* seqCollection, int kmerSize)
{
	Timer timer("PopBubbles");
	printf("removing loops\n");
	int numPopped = 0;

	// Set the cutoffs
	const unsigned int expectedBubbleSize = 2*(kmerSize + 1);
	const unsigned int maxNumBranches = 3;
		
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		// Skip sequences that have already been deleted	
		if(iter->isFlagSet(SF_DELETE))
		{		
			continue;
		}

		// Get the extensions for this sequence, this function populates the extRecord structure
		ExtensionRecord extRec;
		bool success = seqCollection->getExtensions(*iter, extRec);
		assert(success);
		
		// Check for ambiguity
		for(int i = 0; i <= 1; ++i)
		{	
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			
			if(extRec.dir[dir].IsAmbiguous())
			{
				// Found a potential bubble, examine each branch
				bool stop = false;
				
				// Create the branch group
				BranchGroup branchGroup(0, dir, maxNumBranches);
				initiateBranchGroup(branchGroup, *iter, extRec.dir[dir], expectedBubbleSize);
				
				// Iterate over the branches
				while(!stop)
				{
					size_t numBranches = branchGroup.getNumBranches();
					
					for(unsigned int j = 0; j < numBranches; ++j)
					{						
						// Get the extensions of this branch
						ExtensionRecord extRec;
						bool success = seqCollection->getExtensions(branchGroup.getBranch(j).getLastSeq(), extRec);
						assert(success);
						
						processBranchGroupExtension(branchGroup, j, branchGroup.getBranch(j).getLastSeq(), extRec);
					}
					
					// At this point all branches should have the same length or one will be a noext
					branchGroup.updateStatus();
					
					// All branches have been extended one sequence, check the stop conditions
					BranchGroupStatus status = branchGroup.getStatus();
					
					// Check if a stop condition was met
					if(status == BGS_TOOLONG || status == BGS_LOOPFOUND || status == BGS_TOOMANYBRANCHES || status == BGS_NOEXT)
					{
						stop = true;
					}
					else if(status == BGS_JOINED)
					{
						collapseJoinedBranches(seqCollection, branchGroup);
						numPopped++;
						stop = true;
					}
					else
					{										
						// the branch is still active, continue
						assert(status == BGS_ACTIVE);
					}
				}
			}
		}
		seqCollection->pumpNetwork();
	}
	printf("Removed %d bubbles\n", numPopped);
	return numPopped;
}

//
// Populate a branch group with the inital branches from a sequence
//
void initiateBranchGroup(BranchGroup& group, const PackedSeq& seq, const SeqExt& extension, size_t maxBubbleSize)
{
	// Generate the vector of sequences that make up this branch
	PSequenceVector extSeqs;
	generateSequencesFromExtension(seq, group.getDirection(), extension, extSeqs);
	assert(extSeqs.size() > 1);
	uint64_t id = 0;
	for(PSequenceVector::iterator seqIter = extSeqs.begin(); seqIter != extSeqs.end(); ++seqIter)
	{
		// Create a new branch and add it to the group
		BranchRecord newBranch(group.getDirection(), maxBubbleSize);
		BranchRecord& addedBranch = group.addBranch(id, newBranch);
		
		// Add the sequence to the branch
		addedBranch.addSequence(*seqIter);
		
		id++;
	}	
}

//
// process an a branch group extension
//
bool processBranchGroupExtension(BranchGroup& group, size_t branchIndex, const PackedSeq& seq, ExtensionRecord extensions)
{
	// Generate the extensions of the branch
	PSequenceVector branchExtSeqs;
	extDirection dir = group.getDirection();
	generateSequencesFromExtension(seq, dir, extensions.dir[dir], branchExtSeqs);
	
	if(branchExtSeqs.size() == 1)
	{
		// single extension
		group.getBranch(branchIndex).addSequence(branchExtSeqs.front());
							
	}
	else if(branchExtSeqs.size() > 1)
	{
		// Start a new branch for the sequences [1..n]
		PSequenceVector::iterator seqIter = branchExtSeqs.begin() + 1;
		
		uint64_t startID = group.getNumBranches();
		for(; seqIter != branchExtSeqs.end(); ++seqIter)
		{
			//printf("adding multip %s\n",bubHR.getHit(k).seq.decode().c_str());
			
			// Start a new branch which is a duplicate of the current branch up to this point
			BranchRecord newBranch(group.getBranch(branchIndex));
			BranchRecord& addedBranch = group.addBranch(startID, newBranch);
			addedBranch.addSequence(*seqIter);
		}
		
		// Add the first sequence (index 0) to the current branch
		group.getBranch(branchIndex).addSequence(branchExtSeqs.front());
	}
	else
	{
	
		// this branch could not be extended, set a flag
		group.setNoExtension();
	}
	
	// Return whether the group is extendable
	return group.isExtendable();
}

//
// Collapse joined paths into a single path
//
void collapseJoinedBranches(ISequenceCollection* seqCollection, BranchGroup& group)
{
	// a join was found, select a branch to keep and remove the rest
	size_t selectedIndex = group.selectBranchToKeep();
	
	size_t numBranches = group.getNumBranches();
	BranchRecord& refRecord = group.getBranch(selectedIndex);
	
	for(size_t i = 0; i < numBranches; ++i)
	{
		// Skip the branch that was selected to keep
		if(i == selectedIndex)
		{
			continue;
		}
		
		BranchRecord& currBranch = group.getBranch(i);
		for(BranchDataIter branchIter = currBranch.getStartIter(); branchIter != currBranch.getEndIter(); ++branchIter)
		{
			// if the sequence is not in the reference branch, remove it
			if(!refRecord.exists(*branchIter))
			{
				//printf("Deleting %s from branch %x\n", branchIter->decode().c_str(), i);
				removeSequenceAndExtensions(seqCollection, *branchIter);	
			}
		}
	}
}
//
// Remove a read and update the extension records of the sequences that extend to it
//
void removeSequenceAndExtensions(ISequenceCollection* seqCollection, const PackedSeq& seq)
{
	// This removes the reverse complement as well
	seqCollection->remove(seq);
	removeExtensionsToSequence(seqCollection, seq, SENSE);
	removeExtensionsToSequence(seqCollection, seq, ANTISENSE);
}

//
//
//
void removeExtensionsToSequence(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection dir)
{
	extDirection oppDir = oppositeDirection(dir);
	
	PackedSeq testSeq(seq);
	
	char extBase = testSeq.rotate(dir, BASES[0]);
	
	for(int i = 0; i < NUM_BASES; i++)
	{	
		char currBase = BASES[i];
		// don't bother checking if the extension exists, just remove it
		// if the sequence didnt have the extension, the operation will do nothing
		testSeq.setLastBase(dir, currBase);
		
		// remove the extension, this removes the reverse complement as well
		seqCollection->removeExtension(testSeq, oppDir, extBase);
	}	
}

//
// Trimming driver function
//
void performTrim(ISequenceCollection* seqCollection)
{
	if (opt::trimLen == 0)
		return;

	int start = 2;
	while(start <= opt::trimLen)
	{
		trimSequences(seqCollection, start);
		start <<= 1;
	}
	
	bool stop = false;
	while(!stop)
	{
		int numRemoved = trimSequences(seqCollection, opt::trimLen);
		if(numRemoved <= 0)
		{
			stop = true;
		}
	}
}

//
// Check the adjacency of a sequence
//
SeqContiguity checkSeqContiguity(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection& outDir)
{
	if(seqCollection->checkFlag(seq, SF_DELETE) || seqCollection->checkFlag(seq, SF_SEEN) )
	{
		return SC_INVALID;
	}
			
	bool child = seqCollection->hasChild(seq);
	bool parent = seqCollection->hasParent(seq);
	
	if(!child && !parent)
	{
		//this sequence is completely isolated
		return SC_ISLAND;
	}
	else if(!child)
	{
		outDir = ANTISENSE;
		return SC_ENDPOINT;
	}
	else if(!parent)
	{
		outDir = SENSE;
		return SC_ENDPOINT;
	}
	else
	{
		// sequence is contiguous
		return SC_CONTIGUOUS;
	}
}

//
// Trimming (error removal) function
//
int trimSequences(ISequenceCollection* seqCollection, int maxBranchCull)
{
	Timer timer("TrimSequences");
	printf("trimming max branch: %d\n", maxBranchCull);	
	int numBranchesRemoved = 0;

	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		
		extDirection dir;
		// dir will be set to the trimming direction if the sequence can be trimmed
		SeqContiguity status = checkSeqContiguity(seqCollection, *iter, dir);

		if(status == SC_INVALID || status == SC_CONTIGUOUS)
		{
			continue;
		}
		else if(status == SC_ISLAND)
		{
			// remove this sequence, it has no extensions
			removeSequenceAndExtensions(seqCollection, *iter);
		}
		// Sequence is trimmable, continue

		// This is a dead-end branch, check it for removal
		BranchRecord currBranch(dir, maxBranchCull);
					
		PackedSeq currSeq = *iter;
		
		while(currBranch.isActive())
		{		
			// Get the extensions for this sequence, this function populates the extRecord structure
			ExtensionRecord extRec;
			bool success = seqCollection->getExtensions(currSeq, extRec);
			assert(success);
			
			// process the extension record and extend the current branch, this function updates currSeq on successful extension
			processLinearExtensionForBranch(currBranch, currSeq, extRec);
		}
		
		// The branch has ended check it for removal, returns true if it was removed
		if(processTerminatedBranchTrim(seqCollection, currBranch))
		{
			numBranchesRemoved++;
		}
		seqCollection->pumpNetwork();
	}
	
	printf("num branches removed: %d\n", numBranchesRemoved);
	return numBranchesRemoved;
}

//
// Process the extension for this branch for the trimming algorithm
// CurrSeq is the current sequence being inspected (the next member to be added to the branch). The extension record is the extensions of that sequence
// After processing currSeq is unchanged if the branch is no longer active or else it is the generated extension
//
bool processLinearExtensionForBranch(BranchRecord& branch, PackedSeq& currSeq, ExtensionRecord extensions)
{
	extDirection dir = branch.getDirection();
	extDirection oppDir = oppositeDirection(dir);
	
	if(branch.doLengthCheck() && ((int)branch.getLength() > (int)branch.getMaxLength())) // Check if the branch has extended past the max trim length
	{
		// Stop the branch
		branch.terminate(BS_TOO_LONG);
	}
	else if(branch.hasLoop())
	{
		branch.terminate(BS_LOOP);
	}
	else if(extensions.dir[oppDir].IsAmbiguous()) // Does this sequence split TO the former node?
	{
		//printf("stopped because of reverse branch\n");
		// There is a reverse ambiguity to this branch, stop the branch without adding the current sequence to it
		branch.terminate(BS_AMBI_OPP);
	}
	else if(!extensions.dir[dir].HasExtension()) 
	{
		// no extenstion, add the current sequence and terminate the branch
		branch.addSequence(currSeq);
		branch.terminate(BS_NOEXT);
	}
	else if(extensions.dir[dir].IsAmbiguous())
	{
		// this branch has an ambiguous extension, add the current sequence and terminate
		branch.addSequence(currSeq);
		branch.terminate(BS_AMBI_SAME);
	}
	else
	{
		// Add the sequence to the branch
		branch.addSequence(currSeq);
		
		// generate the new current sequence from the extension
		//printf("currseq: %s ", currSeq.decode().c_str());
		PSequenceVector newSeqs;

		generateSequencesFromExtension(currSeq, dir, extensions.dir[dir], newSeqs);
		assert(newSeqs.size() == 1);
		currSeq = newSeqs.front();
		//printf("newseq: %s \n", currSeq.decode().c_str());
	}
	
	return branch.isActive();
}

//
//
//
bool processTerminatedBranchTrim(ISequenceCollection* seqCollection, BranchRecord& branch)
{
	assert(!branch.isActive());
	//printf("	branch has size: %d\n", branchElements.size());
	if(branch.getLength() > 0 && branch.getState() != BS_TOO_LONG)
	{		
		BranchDataIter endIter  = branch.getEndIter();
		for(BranchDataIter bIter = branch.getStartIter(); bIter != endIter; bIter++)
		{
			removeSequenceAndExtensions(seqCollection, *bIter);
		}
		return true;
	}	
	else
	{
		return false;
	}
}

//
// Assembly function
//
void assemble(ISequenceCollection* seqCollection, int /*readLen*/, int /*kmerSize*/, IFileWriter* fileWriter)
{
	Timer timer("Assemble");
	
	int numAssembled = 0;
	int contigID = 0;
	
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		
		extDirection dir;
		// dir will be set to the trimming direction if the sequence can be trimmed
		SeqContiguity status = checkSeqContiguity(seqCollection, *iter, dir);

		if(status == SC_INVALID || status == SC_CONTIGUOUS)
		{
			continue;
		}
		else if(status == SC_ISLAND)
		{
			// this should never happen, island sequences would have been removed by this point
			assert(false);
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
			processLinearExtensionForBranch(currBranch, currSeq, extRec);
		}
		
		// The branch has ended, build the contig
		Sequence contig;
		processTerminatedBranchAssemble(seqCollection, currBranch, contig);
		
		// Output the contig
		fileWriter->WriteSequence(contig, contigID, 0);
		
		contigID++;				
		numAssembled++;

		seqCollection->pumpNetwork();
	}
	
	printf("num branches assembled: %d\n", numAssembled);
}

//
//
//
void processTerminatedBranchAssemble(ISequenceCollection* seqCollection, BranchRecord& branch, Sequence& outseq)
{
	assert(!branch.isActive());
	//printf("	branch has size: %d\n", branchElements.size());
	
	// the only acceptable condition for the termination of an assembly is a noext or a loop
	assert(branch.getState() == BS_NOEXT || branch.getState() == BS_LOOP);
	
	// Mark all the sequences in the branch as seen
	BranchDataIter endIter  = branch.getEndIter();
	for(BranchDataIter bIter = branch.getStartIter(); bIter != endIter; bIter++)
	{
		seqCollection->setFlag(*bIter, SF_SEEN);
	}
	
	// Assemble the contig
	branch.buildContig(outseq);
}

// Write the sequences out to a file
//
void outputSequences(const char* filename, ISequenceCollection* pSS)
{
	FastaWriter writer(filename);
	SequenceCollectionIterator endIter  = pSS->getEndIter();
	int64_t count = 0;
	for(SequenceCollectionIterator iter = pSS->getStartIter(); iter != endIter; ++iter)
	{
		if(!pSS->checkFlag(*iter, SF_DELETE))
		{
			writer.WriteSequence(iter->decode(), count, 0.0f);
			count++;
		}
	}	
}

};
