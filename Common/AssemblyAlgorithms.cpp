#include <iostream>
#include "AssemblyAlgorithms.h"
#include "Options.h"
#include "SequenceCollectionHash.h"


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
}

//
// Generate the adjacency information for each sequence in the collection
//
void generateAdjacency(ISequenceCollection* seqCollection)
{
	
	printf("generating adjacency info\n");
	double start = std::clock();

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
	
	double ticks = std::clock() - start;
	double time = (ticks) / (double)CLOCKS_PER_SEC;		
	printf("adjacency set %d bases in %lf seconds\n", numBasesSet, time);
}

//
//
//
void splitAmbiguous(ISequenceCollection* seqCollection)
{
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

void popBubbles(ISequenceCollection* seqCollection, int kmerSize)
{
	printf("removing loops\n");
	int count = 0;
	int numPopped = 0;

	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
	
		if(iter->isFlagSet(SF_DELETE))
		{		
			continue;
		}
				
		if(count % 100000 == 0)
		{
			printf("checked %d for bubbles\n", count);
		}
		count++;
		
		// Set the cutoffs
		unsigned int expectedBubbleSize = 2*(kmerSize + 1);
		const unsigned int maxNumBranches = 3;

		for(int i = 0; i <= 1; i++)
		{
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			HitRecord initHR = calculateExtension(seqCollection, *iter, dir);

			if(initHR.getNumHits() > 1)
			{
				// Found a potential bubble, examine each branch
				bool stop = false;
						
				std::vector<Branch> branches;
				
				// Add the initial sequences to the branch
				for(int j = 0; j < initHR.getNumHits(); ++j)
				{
					Branch newBranch;
					//printf("Adding start: %s\n", initHR.getHit(j).seq.decode().c_str());
					
					newBranch.AddSequence(initHR.getHit(j).seq);
					branches.push_back(newBranch);
				}

				
				bool bubble = false;
				// Iterate over the branches
				while(!stop)
				{
					int numBranches = branches.size();
					for(int j = 0; j < numBranches; ++j)
					{												
						// Check the extension of this sequence
						HitRecord bubHR = calculateExtension(seqCollection, branches[j].lastSeq, dir);
						//printf("currSeq: %s size: %d num: %d\n", branches[j].lastSeq.decode().c_str(), branches[j].seqSet.size(), branches.size());
						
						bool uniqueEntry = false;
						if(bubHR.getNumHits() == 1)
						{
							// single extension
							//printf("adding single %s\n", bubHR.getFirstHit().seq.decode().c_str());
							uniqueEntry = branches[j].AddSequence(bubHR.getFirstHit().seq);
												
						}
						else if(bubHR.getNumHits() > 1)
						{
							// Start a new branch for the remaining sequences
							for(int k = 1; k < bubHR.getNumHits(); k++)
							{
								//printf("adding multip %s\n",bubHR.getHit(k).seq.decode().c_str());
								// Start a new branch which is a duplicate of the current branch up to this point
								
								// Always unique so we don't have to check the return
								Branch newBranch(branches[j]);
								newBranch.AddSequence(bubHR.getHit(k).seq);							
								branches.push_back(newBranch);	
							}
										
							// Add the first base to the current branch
							//printf("adding multip %s\n",bubHR.getFirstHit().seq.decode().c_str());
							uniqueEntry = branches[j].AddSequence(bubHR.getFirstHit().seq);
						}
						
						// If there was no extension or we tried to add a duplicate key, return
						if(bubHR.getNumHits() == 0 || !uniqueEntry)
						{
							// no ext, terminate
							stop = true;
							bubble = false;
							break;
						}
					}
					
					// Check the stop conditions	

					
					// Check if all the branches have joined
					bool joinFound = false;
					for(unsigned int k = 0; k < branches.size(); k++)
					{
						bool allSame = true;
						PackedSeq testSeq = branches[k].lastSeq;
						for(unsigned int l = 0; l < branches.size(); l++)
						{
							// Skip the branch that is the same as the current
							if(k == l)
							{
								continue;
							}
							
							if(branches[l].seqSet.find(testSeq) == branches[l].seqSet.end())
							{
								// Sequence not found
								allSame = false;
								break;
							}	
						}
						
						if(allSame)
						{
							joinFound = true;
							break;
						}
					}

					if(joinFound)
					{
						// The paths have merged back together, this is a bubble
						stop = true;
						bubble = true;
					}
					
					// Check stop conditions
					if(branches.size() > maxNumBranches)
					{
						stop = true;
						bubble = false;
					}
					
					for(unsigned int k = 0; k < branches.size(); k++)
					{
						if(branches[k].seqSet.size() > expectedBubbleSize)
						{
							stop = true;
							bubble = false;
						}
					}	
				}
				/*
				printf("Bubble found of size %d\n", branches[0].size());
				for(unsigned int k = 0; k < branches.size(); k++)
				{
					printf("Path %d:\n", k);
					for(PSequenceVectorIterator bubIter = branches[k].begin(); bubIter != branches[k].end(); bubIter++)
					{
						printf("	%s\n", bubIter->decode().c_str());
					}
				}				
				*/
				// Was a bubble found?
				if(bubble)
				{
					//printf("BRANCHES %d %d\n", branches.size(), branches.front().size());
					// Remove the bubble (arbitrary decision on which one for now)
					//printf("pop!\n");
					for(unsigned int k = 1; k < branches.size(); k++)
					{
						int removeIndex = k;
						
						// Stop removal at the second last sequence
						for(PSeqSet::iterator bubIter = branches[removeIndex].seqSet.begin(); bubIter != branches[removeIndex].seqSet.end(); bubIter++)
						{
							if(*bubIter != branches[removeIndex].lastSeq && branches[0].seqSet.find(*bubIter) == branches[0].seqSet.end())
							{
								//printf("Deleting (%zu): %s\n", branches[removeIndex].seqSet.size(), bubIter->decode().c_str());
								removeSequenceAndExtensions(seqCollection, *bubIter);
							}
						}
					}
					numPopped++;

					//printf("Popped %zu\n", branches[0].seqSet.size());
				}
			}
		}
		seqCollection->pumpNetwork();
	}
	printf("Removed %d bubbles\n", numPopped);	
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
		
	for(int i = 0; i < NUM_BASES; i++)
	{	
		char currBase = BASES[i];
		// does this sequence have an extension to the deleted seq?
		ResultPair hasExt  = seqCollection->checkExtension(seq, dir, currBase);
		if(hasExt.forward || hasExt.reverse)
		{
			PackedSeq tempSeq(seq);	
			// generate the sequence that the extension is to
			char extBase = tempSeq.rotate(dir, currBase);				
			// remove the extension, this removes the reverse complement as well
			seqCollection->removeExtension(tempSeq, oppDir, extBase);
		}
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
// Check whether a sequence can be trimmed
//
TrimStatus checkSeqForTrim(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection& outDir)
{
	if(seqCollection->checkFlag(seq, SF_DELETE))
	{
		return TS_NOTRIM;
	}
			
	bool child = seqCollection->hasChild(seq);
	bool parent = seqCollection->hasParent(seq);
	
	if(!child && !parent)
	{
		//this sequence is completely isolated
		return TS_ISLAND;
	}
	else if(!child)
	{
		outDir = ANTISENSE;
		return TS_TRIMMABLE;
	}
	else if(!parent)
	{
		outDir = SENSE;
		return TS_TRIMMABLE;
	}
	else
	{
		// sequence is contiguous
		return TS_NOTRIM;	
	}
}

//
// Trimming (error removal) function
//
int trimSequences(ISequenceCollection* seqCollection, int maxBranchCull)
{
	printf("trimming max branch: %d\n", maxBranchCull);	
	int numBranchesRemoved = 0;

	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		
		extDirection dir;
		// dir will be set to the trimming direction if the sequence can be trimmed
		TrimStatus status = checkSeqForTrim(seqCollection, *iter, dir);

		if(status == TS_NOTRIM)
		{
			continue;
		}
		else if(status == TS_ISLAND)
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
			
			// this call should always succeed
			if(!success)
			{
				printf("seq %s was not found!\n", currSeq.decode().c_str());
			}
			assert(success);
			
			// process the extension record and extend the current branch, this function updates currSeq on successful extension
			processExtensionForBranchTrim(currBranch, currSeq, extRec);
		}
		
		// The branch has ended check it for removal, returns true if it was removed
		if(processTerminatedBranch(seqCollection, currBranch))
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
bool processExtensionForBranchTrim(BranchRecord& branch, PackedSeq& currSeq, ExtensionRecord extensions)
{
	extDirection dir = branch.getDirection();
	extDirection oppDir = oppositeDirection(dir);
	
	if(branch.getLength() > branch.getMaxLength()) // Check if the branch has extended past the max trim length
	{
		// Stop the branch
		branch.terminate(BS_TOO_LONG);
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
bool processTerminatedBranch(ISequenceCollection* seqCollection, BranchRecord& branch)
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
void assemble(ISequenceCollection* seqCollection, int readLen, int kmerSize, IFileWriter* fileWriter)
{

	int noext = 0;
	int ambiext = 0;

	printf("starting assembly\n");
	int count = 0;
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	
	int contigID = 0;
	

	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		if(!seqCollection->checkFlag(*iter, SF_SEEN) && !seqCollection->checkFlag(*iter, SF_DELETE))
		{			
			// Is this sequence a branch endpoint?
			bool doAssembly = false;
			extDirection dir = SENSE;
			if(!seqCollection->hasParent(*iter))
			{
				doAssembly = true;
				dir = SENSE;
			}
			else if(!seqCollection->hasChild(*iter))
			{
				doAssembly = true;
				dir = ANTISENSE;
			}

			int sumMultiplicity = 0;
			if(doAssembly)
			{
				// the record of extensions			
				HitVector extensions[2];
									
				//for(int i = 0; i <= 1; i++)
				{
					bool stop = false;
					PackedSeq currSeq = *iter;
					SeqRecord loopCheck;			
					
					while(!stop)
					{
						//assert(!seqCollection->checkFlag(currSeq, SF_SEEN));
						// Mark the flag for the selected sequence
						sumMultiplicity += seqCollection->getMultiplicity(currSeq);
						seqCollection->setFlag(currSeq, SF_SEEN);
										
						HitRecord hr = calculateExtension(seqCollection, currSeq, dir);
						if(hr.getNumHits() == 0)
						{
							// no ext
							stop = true;
							noext++;
						}
						else if(hr.getNumHits() == 1)
						{
							// good ext
							currSeq = hr.getFirstHit().seq;
							
							if(loopCheck.contains(currSeq))
							{
								stop = true;
							}
							else
							{
								//printf("good ext (%s)\n", currSeq.decode().c_str());
								extensions[dir].push_back(hr.getFirstHit());
							}
						}
						else
						{
							// ambi ext
							stop = true;
							ambiext++;
						}
						
						loopCheck.addSequence(currSeq);
					}				
				}
				Sequence contig = BuildContig(seqCollection, extensions, *iter, contigID, readLen, kmerSize);
				
				// is this contig worth outputting?
				//if(contig.length() >= 100)
				{
					count++;
					//sprintf(buffer, ">%d\n%s\n", count, contig.c_str());
					//printf("%s", buffer);
					//MPI_Status status;
					//MPI_File_write(handle, buffer, numChars, MPI::CHAR, &status);
					double numReads = (double)sumMultiplicity / (readLen - kmerSize + 1);
					double coverage = (numReads * (double)readLen) /  contig.length();
					fileWriter->WriteSequence(contig, contigID, coverage);
					contigID++;
				}					
			}
		}
		seqCollection->pumpNetwork();
	}
	//MPI_File_close(&handle);
	printf("noext: %d, ambi: %d\n", noext, ambiext);	
	
}

//
// Build the contig from the extension information
//
Sequence BuildContig(ISequenceCollection* /*seqCollection*/,
		HitVector* extensions, const PackedSeq& originalSeq,
		int /*contigID*/, int /*readLen*/, int /*kmerSize*/)
{
	Sequence contig;
	
	contig.reserve(originalSeq.getSequenceLength() + extensions[0].size() + extensions[1].size());
	
	int position = 0;
	// output the contig
	// output all the antisense extensions
	for(HitVector::reverse_iterator asIter = extensions[1].rbegin(); asIter != extensions[1].rend(); asIter++)
	{
		contig.append(1, asIter->seq.getFirstBase());
		position++;
	}
	
	// output the current sequence itself
	contig.append(originalSeq.decode());
	position++;
	
	// output the sense extensions
	for(HitVector::iterator sIter = extensions[0].begin(); sIter != extensions[0].end(); sIter++)
	{
		contig.append(1, sIter->seq.getLastBase());
		position++;		
	}	
	return contig;
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
