#include "AssemblyAlgorithms.h"

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
	
	while(reader->isGood())
	{
		PackedSeq seq = reader->ReadSequence();
		
		assert(kmerSize <= seq.getSequenceLength());
		//assert(seq.getSequenceLength() == pairSeq.getSequenceLength());
		
		//int l = pairSeq.getSequenceLength();
		
		for(int i = 0; i < seq.getSequenceLength() - kmerSize  + 1; i++)
		{
			PackedSeq sub = seq.subseq(i, kmerSize);
			SeqID currID(idNum, i);
			
			sub.addID(currID);
			
			//PackedSeq pairSub = pairSeq.subseq(l - kmerSize - i, kmerSize);

			// Add the sequence to the sequence collection
			seqCollection->add(sub);

			
			if(count % 1000000 == 0)
			{
				printf("read %d sequences\n", count);
			}
			count++;
		}
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
	int count = 0;
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
			for(int j = 0; j < NUM_BASES; j++)
			{
				char currBase = BASES[j];
				PackedSeq testSeq(*iter);
				testSeq.rotate(dir, currBase);
				
				if(seqCollection->exists(testSeq))
				{
					extension.SetBase(currBase);
				}
			}
			seqCollection->setExtension(*iter, dir, extension);			
		}
		
		
		//iter->printExtension();
		seqCollection->pumpNetwork();
	}
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
		/*
		PackedSeq testSeq("GCTGGCACCACCACCCCTGGCCACCCCAG");
		if(*iter != testSeq && *iter != reverseComplement(testSeq))
		{
			continue;
		}
		*/
		
		if(iter->isFlagSet(SF_DELETE))
		{
			continue;
		}
				
		if(count % 1000000 == 0)
		{
			printf("checked %d for bubbles\n", count);
		}
		count++;
		
		// Set the cutoffs
		unsigned int expectedBubbleSize = 2*(kmerSize + 1);
		const unsigned int maxNumBranches = 4;
		
		for(int i = 0; i <= 1; i++)
		{
			extDirection dir = (i == 0) ? SENSE : ANTISENSE;
			HitRecord initHR = calculateExtension(seqCollection, *iter, dir);

			if(initHR.getNumHits() > 1)
			{
				//printf("Found potential bubble\n");
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
						//printf("Curr seq: %s\n", currSeq.decode().c_str());
												
						// Check the extension of this sequence
						HitRecord bubHR = calculateExtension(seqCollection, branches[j].lastSeq, dir);
						if(bubHR.getNumHits() == 1)
						{
							// single extension
							//printf("adding single %s\n", bubHR.getFirstHit().seq.decode().c_str());
							branches[j].AddSequence(bubHR.getFirstHit().seq);
												
						}
						else if(bubHR.getNumHits() > 1)
						{
							// Start a new branch for the remaining sequences
							for(int k = 1; k < bubHR.getNumHits(); k++)
							{
								//printf("adding multip %s\n",bubHR.getHit(k).seq.decode().c_str());
								// Start a new branch which is a duplicate of the current branch up to this point
								Branch newBranch(branches[j]);
								newBranch.AddSequence(bubHR.getHit(k).seq);							
								branches.push_back(newBranch);	
							}
										
							// Add the first base to the current branch
							//printf("adding multip %s\n",bubHR.getFirstHit().seq.decode().c_str());
							branches[j].AddSequence(bubHR.getFirstHit().seq);
						}
						else
						{
							// no ext, terminate
							stop = true;
							bubble = false;
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
					
					if(branches[0].seqSet.size() > expectedBubbleSize || branches.size() > maxNumBranches)
					{
						//printf("bubble too long stop\n");
						stop = true;
						bubble = false;
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
					for(unsigned int k = 1; k < branches.size(); k++)
					{
						int removeIndex = k;
						
						// Stop removal at the second last sequence
						for(PSeqSet::iterator bubIter = branches[removeIndex].seqSet.begin(); bubIter != branches[removeIndex].seqSet.end(); bubIter++)
						{
							if(*bubIter != branches[removeIndex].lastSeq)
							{
								printf("Deleting (%lu): %s\n", branches[removeIndex].seqSet.size(), bubIter->decode().c_str());
								removeSequenceAndExtensions(seqCollection, *bubIter);
							}
						}
					}
					numPopped++;
					printf("Popped %lu\n", branches[0].seqSet.size());
				}

		
			}
		}
		seqCollection->pumpNetwork();
	}
	printf("Removed %d bubbles\n", numPopped);	
}

//
//
//
HitRecord calculateExtension(const ISequenceCollection* seqCollection,
		const PackedSeq& currSeq, extDirection dir)
{
	
	// Create the return structure
	HitRecord hitRecord;
	
	// Check for the existance of the 4 possible extensions
	for(int i  = 0; i < NUM_BASES; i++)
	{
		char currBase = BASES[i];
		
		ResultPair hasExt = seqCollection->checkExtension(currSeq, dir, currBase);
			
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
void performTrim(ISequenceCollection* seqCollection, int readLen, int kmerSize)
{
	int start = 2;
	int maxBranch =  6 * (readLen - kmerSize + 1);
	
	
	while(start <= maxBranch)
	{
		trimSequences(seqCollection, start);
		start <<= 1;
	}
	
	bool stop = false;
	while(!stop)
	{
		int numRemoved = trimSequences(seqCollection, maxBranch);
		if(numRemoved <= 0)
		{
			stop = true;
		}
	}
	
	/*	

	// Now trim at the max branch length
	for(int i = 0; i < 2; i++)
	{
		trimSequences(seqCollection, maxBranch);
	}
	
	// finally, trim at the min contig length
	bool stop = false;
	while(!stop)
	{
		int numRemoved = trimSequences(seqCollection, 100);
		if(numRemoved <= 0)
		{
			stop = true;
		}
	}
	*/
}

//
// Trimming (error removal) function
//

int trimSequences(ISequenceCollection* seqCollection, int maxBranchCull)
{
	const int MAX_DEAD_LENGTH = maxBranchCull;
	printf("trimming max branch: %d\n", maxBranchCull);	
	int numBranchesRemoved = 0;
	int count = 0;
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		if(count % 1000000 == 0)
		{
			printf("trimmed: %d\n", count);
		}
		count++;

		if(iter->isFlagSet(SF_DELETE))
		{
			continue;
		}
				
		bool child = seqCollection->hasChild(*iter);
		bool parent = seqCollection->hasParent(*iter);
		
		extDirection dir;
		if(!child && !parent)
		{
			// remove this sequence, it has no extensions
			removeSequenceAndExtensions(seqCollection, *iter);
			continue;
		}
		else if(!child)
		{
			dir = ANTISENSE;
		}
		else if(!parent)
		{
			dir = SENSE;
		}
		else
		{
			continue;	
		}

		PSequenceVector branchElements;
		
		extDirection oppositeDir = oppositeDirection(dir);
					
		PackedSeq currSeq = *iter;
		SeqRecord loopCheck;
		bool stop = false;
		
		while(!stop)
		{		
				
			HitRecord hr = calculateExtension(seqCollection, currSeq, dir);
			HitRecord oppHr = calculateExtension(seqCollection, currSeq, oppositeDir);
			
			if((int)branchElements.size() == MAX_DEAD_LENGTH + 1 )
			{
				// no ext
				stop = true;
				//printf("stopped because of too long: %d\n", branchElements.size());
			}
			else if(oppHr.getNumHits() > 1)
			{
				//printf("stopped because of reverse branch\n");
				stop = true;	
			}
			else if(hr.getNumHits() == 0 || hr.getNumHits() > 1)
			{
				branchElements.push_back(currSeq);
				//printf("stopped because of noext/ambi branch\n");
				stop = true;
			}
			else
			{
				branchElements.push_back(currSeq);
					
				// good ext
				currSeq = hr.getFirstHit().seq;
			}
		}
		
		//printf("	branch has size: %d\n", branchElements.size());
		if((int)branchElements.size() <= MAX_DEAD_LENGTH && (int)branchElements.size() > 0)
		{
			//printf("		removing\n");
			numBranchesRemoved++;
			for(PSequenceVectorIterator bIter = branchElements.begin(); bIter != branchElements.end(); bIter++)
			{
				removeSequenceAndExtensions(seqCollection, *bIter);
			}
		}
		
		seqCollection->pumpNetwork();
	}
	
	
	
	printf("seqs after trimming: %d\n", seqCollection->count());
	printf("num branches removed: %d\n", numBranchesRemoved);
	return numBranchesRemoved;
}


bool doPrint = false;
//
// Assembly function
//
void assemble(ISequenceCollection* seqCollection, int readLen, int kmerSize)
{
	// create file writer
	FastaWriter writer("contigs.fa");


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

			
			// There are 2 cases in which we should extend a read:
			// 1) It is at the endpoint of a branch (either no parent or child extension)
			// 2) It is the first read past a branch
			
			// Is this sequence a branch endpoint?
			bool doAssembly = false;
			extDirection dir;
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
			else
			{
				/*
				// Not a branch endpoint, check if it is the first node after an ambiguous node
				HitRecord senseHR = calculateExtension(seqCollection, *iter, SENSE);
				if(senseHR.getNumHits() == 1)
				{
					// the next node is umambiguous
					HitRecord inHR = calculateExtension(seqCollection, senseHR.getFirstHit().seq, ANTISENSE);
					if(inHR.getNumHits() > 1)
					{
						doAssembly = true;
					}
				}
				
				// Get the next 
				HitRecord antisenseHR = calculateExtension(seqCollection, *iter, ANTISENSE);
				if(antisenseHR.getNumHits() == 1)
				{
					// the next node is umambiguous
					HitRecord inHR = calculateExtension(seqCollection, antisenseHR.getFirstHit().seq, SENSE);
					if(inHR.getNumHits() > 1)
					{
						doAssembly = true;
					}
				}
				*/
			}

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
					writer.WriteSequence(contig, contigID);
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
	std::set<int64_t> readsAligned;
	
	contig.reserve(originalSeq.getSequenceLength() + extensions[0].size() + extensions[1].size());
	
	int position = 0;
	// output the contig
	// output all the antisense extensions
	for(HitVector::reverse_iterator asIter = extensions[1].rbegin(); asIter != extensions[1].rend(); asIter++)
	{
		contig.append(1, asIter->seq.getFirstBase());
		
		//PrintAlignmentForSeq(seqCollection, readsAligned, asIter->seq, contigID, position, readLen, kmerSize);

		position++;
	}
	
	// output the current sequence itself
	contig.append(originalSeq.decode());
	//PrintAlignmentForSeq(seqCollection, readsAligned, originalSeq, contigID, position, readLen, kmerSize);
	position++;
	
	// output the sense extensions
	for(HitVector::iterator sIter = extensions[0].begin(); sIter != extensions[0].end(); sIter++)
	{
		contig.append(1, sIter->seq.getLastBase());
		//PrintAlignmentForSeq(seqCollection, readsAligned, sIter->seq, contigID, position, readLen, kmerSize);
		position++;		
	}	
	return contig;
}

void PrintAlignmentForSeq(ISequenceCollection* seqCollection, std::set<int64_t>& readsAligned, const PackedSeq& seq, int contigID, int position, int readLen, int kmerSize)
{
	IDList fwdIDs = seqCollection->getIDs(seq);
	PrintAlignment(fwdIDs, readsAligned, contigID, position, seq.decode(), false, readLen, kmerSize);
	
	PackedSeq rcSeq = reverseComplement(seq);
	IDList revIDs = seqCollection->getIDs(rcSeq);
	PrintAlignment(revIDs, readsAligned, contigID, position, rcSeq.decode(), true, readLen, kmerSize);	
}

void PrintAlignment(const IDList& ids, std::set<int64_t>& readsAligned,
		int /*contig*/, int position, const Sequence& /*s*/,
		bool isRC, int readLen, int kmerSize)
{
	for(IDList::const_iterator iter = ids.begin(); iter != ids.end(); iter++)
	{
		int64_t numID = iter->first;
		if(readsAligned.find(numID) == readsAligned.end())
		{
			// Calculate the position of the actual read, not the kmer
			int truePosition;
			if(!isRC)
			{
				truePosition = position - iter->second;
			}
			else
			{
				truePosition = position - (readLen - kmerSize - iter->second);
			}
			assert(false);
			//alignments << iter->first << " " << contig << " " << truePosition << " " << isRC << " " << s << std::endl;
			readsAligned.insert(numID);
		}
	}
}

//
// Write the sequences out to a file
//
void outputSequences(const char* filename, ISequenceCollection* pSS)
{
	FastaWriter writer(filename);
	SequenceCollectionIterator endIter  = pSS->getEndIter();
	int count = 0;
	for(SequenceCollectionIterator iter = pSS->getStartIter(); iter != endIter; ++iter)
	{
		if(!pSS->checkFlag(*iter, SF_DELETE))
		{
			writer.WriteSequence(*iter, count);
			count++;
		}
	}	
}
