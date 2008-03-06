#include "AssemblyAlgorithms.h"

//
// Function to load sequences into the collection
//
void loadSequences(ISequenceCollection* seqCollection, std::string fastaFile, int readLength, int kmerSize)
{
	FastaReader* reader = new FastaReader(fastaFile.c_str());
	
	int count = 0;
	
	// Read the sequences and add them to the network sequence space
	while(reader->isGood())
	{
		PackedSeq seq = reader->ReadSequence();		
		assert(kmerSize <= seq.getSequenceLength());
		
		for(int i = 0; i < seq.getSequenceLength() - kmerSize  + 1; i++)
		{
			PackedSeq sub = seq.subseq(i, kmerSize);
			
			// Add the sequence to the network space
			seqCollection->add(sub);
			
			if(count % 100000 == 0)
			{
				printf("read %d sequences\n", count);
			}
			count++;
		}
		
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
	PSequenceVectorIterator endIter  = seqCollection->getEndIter();
	for(PSequenceVectorIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		if(count % 100000 == 0)
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
HitRecord calculateExtension(ISequenceCollection* seqCollection, const PackedSeq& currSeq, extDirection dir)
{
	
	// Create the return structure
	HitRecord hitRecord;
	
	// Check for the existance of the 4 possible extensions
	for(int i  = 0; i < NUM_BASES; i++)
	{
		char currBase = BASES[i];
		
		// SLOW
		bool hasExt = seqCollection->checkExtension(currSeq, dir, currBase);
			
		// Does this sequence have an extension?
		if(hasExt)
		{
			PackedSeq extSeq(currSeq);
			extSeq.rotate(dir, currBase);
			
			// is there a forward extension?
			if(hasExt)
			{
				hitRecord.addHit(extSeq, false);
			}
			//else
			//{
				// extension is of the reverse complement
				//hitRecord.addHit(extSeq, true);	
			//}
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
	
	// Remove this sequence as an extension to the adjacent sequences
	for(int i = 0; i <= 1; i++)
	{
		extDirection dir = (i == 0) ? SENSE : ANTISENSE;
		extDirection oppDir = oppositeDirection(dir);	
			
		for(int i = 0; i < NUM_BASES; i++)
		{	
			char currBase = BASES[i];
			// does this sequence have an extension to the deleted seq?
			bool hasExt  = seqCollection->checkExtension(seq, dir, currBase);
			if(hasExt)
			{
				PackedSeq tempSeq(seq);	
				// generate the sequence that the extension is to
				char extBase = tempSeq.rotate(dir, currBase);				
				// remove the extension, this removes the reverse complement as well
				seqCollection->removeExtension(tempSeq, oppDir, extBase);
			}
		}
	}
}

//
// Trimming driver function
//
void performTrim(ISequenceCollection* seqCollection, int readLen, int kmerSize)
{
	int start = 2;
	int maxBranch =  4 * (readLen - kmerSize + 1);
	
	
	while(start <= maxBranch)
	{
		trimSequences(seqCollection, start);
		start <<= 1;
	}
	
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
	PSequenceVectorIterator endIter  = seqCollection->getEndIter();
	for(PSequenceVectorIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		
		if(count % 100000 == 0)
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
			
			if(branchElements.size() == MAX_DEAD_LENGTH + 1 )
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
		if(branchElements.size() <= MAX_DEAD_LENGTH && branchElements.size() > 0)
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

//
// Assembly function
//
void assemble(ISequenceCollection* seqCollection)
{
	// create file writer
	FastaWriter writer("contigs.fa");
	int noext = 0;
	int ambiext = 0;

	printf("starting assembly\n");
	int count = 0;
	PSequenceVectorIterator endIter  = seqCollection->getEndIter();
	for(PSequenceVectorIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		if(!seqCollection->checkFlag(*iter, SF_SEEN) && !seqCollection->checkFlag(*iter, SF_DELETE))
		{
			// the record of extensions			
			PSequenceVector extensions[2];
							
			for(int i = 0; i <= 1; i++)
			{
				bool stop = false;
				extDirection dir = (i == 0) ? SENSE : ANTISENSE;
				PackedSeq currSeq = *iter;
				SeqRecord loopCheck;			
				
				while(!stop)
				{
					
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
							extensions[i].push_back(currSeq);
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
			
			Sequence contig = BuildContig(extensions, *iter);
			
			// is this contig worth outputting?
			if(contig.length() >= 100)
			{
				count++;
				//const int bufferSize = 10*1024;
				//char buffer[bufferSize];
				//sprintf(buffer, ">%d\n%s\n", count, contig.c_str());
				//printf("%s", buffer);
				//MPI_Status status;
				//MPI_File_write(handle, buffer, numChars, MPI::CHAR, &status);
				writer.WriteSequence(contig);
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
Sequence BuildContig(PSequenceVector* extensions, const PackedSeq& originalSeq)
{
	Sequence contig;
	contig.reserve(originalSeq.getSequenceLength() + extensions[0].size() + extensions[1].size());
	
	// output the contig
	// output all the antisense extensions
	for(PSequenceVector::reverse_iterator asIter = extensions[1].rbegin(); asIter != extensions[1].rend(); asIter++)
	{
		contig.append(1, asIter->getFirstBase());
	}
	
	// output the current sequence itself
	contig.append(originalSeq.decode());
	
	// output the sense extensions
	for(PSequenceVector::iterator sIter = extensions[0].begin(); sIter != extensions[0].end(); sIter++)
	{
		contig.append(1, sIter->getLastBase());
	}	
	return contig;
}
