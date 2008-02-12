#include <stdio.h>

#include <vector>
#include <stdio.h>
#include "Abyss.h"
#include "CommonUtils.h"

int main(int argv, char** argc)
{	
	if(argv < 4 || argc[1] == "--help")
	{
		printUsage();
		exit(1);
	}
	
	std::string fastaFile = argc[1];
	int kmerSize = atoi(argc[2]);
	int numTrims = atoi(argc[3]);

	// Set up coordinate space
	int maxC = kmerSize - 1;
	Coord4 minCoords = {0, 0, 0, 0};
	Coord4 maxCoords = {maxC, maxC, maxC, maxC};
	
	// Load the phase space
	PhaseSpace* pPS = new PhaseSpace(kmerSize, minCoords, maxCoords);
	
	// Open file for read
	FastaReader* reader = new FastaReader(fastaFile.c_str());
	
	// Load phase space
	int count = 0;

	
	while(reader->isGood())
	{
		PackedSeq seq = reader->ReadSequence();		
		assert(kmerSize <= seq.getSequenceLength());
		
		for(int i = 0; i < seq.getSequenceLength() - kmerSize  + 1; i++)
		{
			PackedSeq sub = seq.subseq(i, kmerSize);
			pPS->addSequence(sub, true);

			count++;
			
			if(count % 1000000 == 0)
			{
				printf("loaded %d sequences\n", count);
			}
		}
	}
	
	// close the reader	
	delete reader;
	
	printf("done sequence load (%d sequences)\n", count);

	// trim the reads
	for(int i = 0; i < numTrims; i++)
	{
		printf("trimming %d/%d\n", i+1, numTrims);
		trimSequences(pPS, minCoords, maxCoords);
		//trimSequences2(pPS, minCoords, maxCoords, i+1);
	}
	
	printf("outputting trimmed reads\n");
	outputSequences("trimmed.fa", pPS, minCoords, maxCoords);
	
	printf("assembling sequences\n");
	assemble(pPS, minCoords, maxCoords);
	
	delete pPS;
	return 0;
}

void outputSequences(const char* filename, PhaseSpace* pPS, Coord4 minCoord, Coord4 maxCoord)
{
	FastaWriter writer(filename);
	for(int x = minCoord.x; x < maxCoord.x; x++) {
		for(int y = minCoord.y; y < maxCoord.y; y++) {
			for(int z = minCoord.z; z < maxCoord.z; z++) {
				for(int w = minCoord.w; w < maxCoord.w; w++)
				{
					Coord4 c = {x,y,z,w};
					
					for(PhaseSpaceBinIter iter = pPS->getStartIter(c); iter != pPS->getEndIter(c); iter++)
					{
						writer.WriteSequence(*iter);
					}	
				}
			}
		}
	}	
}

int noext = 0;
int ambiext = 0;
int count = 0;
	
Sequence BuildContig(PSequenceVector& extensions, PackedSeq& originalSeq, extDirection dir)
{
	Sequence contig;
	contig.reserve(originalSeq.getSequenceLength() + extensions.size());
	
	if(dir == SENSE)
	{
		contig.append(originalSeq.decode());
		for(PSequenceVector::iterator sIter = extensions.begin(); sIter != extensions.end(); sIter++)
		{
			contig.append(1, sIter->getLastBase());
		}
	}
	else
	{
		for(PSequenceVector::reverse_iterator asIter = extensions.rbegin(); asIter != extensions.rend(); asIter++)
		{
			contig.append(1, asIter->getFirstBase());
		}
		
		// output the current sequence itself
		contig.append(originalSeq.decode());		
	}	
	return contig;		
}

void assemble(PhaseSpace* pPS, Coord4 minCoord, Coord4 maxCoord)
{	
	std::list<PackedSeq> seqStarts;
	
	for(int x = minCoord.x; x < maxCoord.x; x++) {
		for(int y = minCoord.y; y < maxCoord.y; y++) {
			for(int z = minCoord.z; z < maxCoord.z; z++) {
				for(int w = minCoord.w; w < maxCoord.w; w++)
				{
					Coord4 c = {x,y,z,w};
					
					for(PhaseSpaceBinIter iter = pPS->getStartIter(c); iter != pPS->getEndIter(c); iter++)
					{
						//if(!pPS->hasChild(*iter) || !pPS->hasParent(*iter))
						{
							seqStarts.push_back(*iter);
						}
					}
				}
			}
		}
	}
	
	printf("num starts: %d\n", seqStarts.size());
	
	// create file writer
	FastaWriter writer("contigs.fa");
	
	while(!seqStarts.empty())
	{
		std::list<PackedSeq>::iterator iter = seqStarts.begin();
	
		if(!pPS->checkSequenceFlag(*iter, SF_SEEN))
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
					pPS->markSequence(currSeq, SF_SEEN);
					pPS->markSequence(reverseComplement(currSeq), SF_SEEN);
									
					HitRecord hr = pPS->calculateExtension(currSeq, dir);
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
						
						/*
						for(int i = 0; i < hr.getNumHits(); i++)
						{
							if(hr.getHit(i).isReverse)
							{
								seqStarts.push_back(reverseComplement(hr.getHit(i).seq));
							}
							else
							{
								seqStarts.push_back(hr.getHit(i).seq);
							}
						}
						*/
					}
					
					loopCheck.addSequence(currSeq);
				}
			}
			
			Sequence contig = BuildContig(extensions, *iter);
			// is this contig worth outputting?
			if(contig.length() >= 100)
			{
				writer.WriteSequence(contig);
			}
		}
		
		// remove this sequence from the todo list
		seqStarts.erase(iter);		
	}
	
	printf("noext: %d, ambi: %d\n", noext, ambiext);
}

void assemble2(PhaseSpace* pPS, Coord4 minCoord, Coord4 maxCoord)
{	
	std::list<branchEnd> seqStarts;
	
	for(int x = minCoord.x; x < maxCoord.x; x++) {
		for(int y = minCoord.y; y < maxCoord.y; y++) {
			for(int z = minCoord.z; z < maxCoord.z; z++) {
				for(int w = minCoord.w; w < maxCoord.w; w++)
				{
					Coord4 c = {x,y,z,w};
					
					for(PhaseSpaceBinIter iter = pPS->getStartIter(c); iter != pPS->getEndIter(c); iter++)
					{
						bool hasChild = pPS->hasChild(*iter);
						bool hasParent = pPS->hasParent(*iter);
						
						if(!hasChild)
						{
							seqStarts.push_back(branchEnd(*iter, ANTISENSE));
						}
						else if(!hasParent)
						{
							seqStarts.push_back(branchEnd(*iter, SENSE));
						}
					}	
				}
			}
		}
	}	
	
	printf("num starts: %d\n", seqStarts.size());
	
	// create file writer
	FastaWriter writer("contigs.fa");
	
	while(!seqStarts.empty())
	{
		std::list<branchEnd>::iterator iter = seqStarts.begin();
		PackedSeq currSeq = iter->first;
		extDirection dir = iter->second;
			
		if(!pPS->checkSequenceFlag(currSeq, SF_SEEN))
		{
			printf("starting from %s\n", currSeq.decode().c_str());
			
			// the record of extensions			
			PSequenceVector extensions;

			bool stop = false;
			SeqRecord loopCheck;			
				
			while(!stop)
			{
					
				// Mark the flag for the selected sequence
				pPS->markSequence(currSeq, SF_SEEN);
				pPS->markSequence(reverseComplement(currSeq), SF_SEEN);
									
				HitRecord hr = pPS->calculateExtension(currSeq, dir);
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
						extensions.push_back(currSeq);
					}
				}
				else
				{
					// ambi ext
					stop = true;
					ambiext++;
					
					printf("exten from %s is ambi\n", currSeq.decode().c_str());
					for(int i = 0; i < hr.getNumHits(); i++)
					{
						if(hr.getHit(i).isReverse)
						{
							seqStarts.push_back(branchEnd(reverseComplement(hr.getHit(i).seq), dir));
						}
						else
						{
							seqStarts.push_back(branchEnd(hr.getHit(i).seq, dir));
						}
					}
					
				}
				
				loopCheck.addSequence(currSeq);
			}
			
			Sequence contig = BuildContig(extensions, iter->first, dir);
			// is this contig worth outputting?
			if(contig.length() >= 100)
			{
				writer.WriteSequence(contig);
			}
		}
		
		// remove this sequence from the todo list
		seqStarts.erase(iter);		
	}
	
	printf("noext: %d, ambi: %d\n", noext, ambiext);
}

Sequence BuildContig(PSequenceVector* extensions, PackedSeq& originalSeq)
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

Sequence assembleSequence(PhaseSpace* pPS, 	PhaseSpaceBinIter sequenceIter)
{		
	Sequence contig;
	return contig;
}

void printUsage()
{
	printf("usage: ABYSS <reads fasta file> <kmer size> <number of trimming steps>\n");	
}



void trimSequences(PhaseSpace* pPS, Coord4 minCoord, Coord4 maxCoord)
{
	for(int x = minCoord.x; x < maxCoord.x; x++) {
		for(int y = minCoord.y; y < maxCoord.y; y++) {
			for(int z = minCoord.z; z < maxCoord.z; z++) {
				for(int w = minCoord.w; w < maxCoord.w; w++)
				{
					Coord4 c = {x,y,z,w};
					
					for(PhaseSpaceBinIter iter = pPS->getStartIter(c); iter != pPS->getEndIter(c); iter++)
					{
						if(!(pPS->hasChild(*iter) && pPS->hasParent(*iter)))
						{
							pPS->removeSequence(*iter);
						}
					}	
				}
			}
		}
	}
}

void trimSequences2(PhaseSpace* pPS, Coord4 minCoord, Coord4 maxCoord, int trimNum)
{	 
	std::vector<branchEnd> deadends;

	for(int x = minCoord.x; x < maxCoord.x; x++) {
		for(int y = minCoord.y; y < maxCoord.y; y++) {
			for(int z = minCoord.z; z < maxCoord.z; z++) {
				for(int w = minCoord.w; w < maxCoord.w; w++)
				{
					Coord4 c = {x,y,z,w};
					
					for(PhaseSpaceBinIter iter = pPS->getStartIter(c); iter != pPS->getEndIter(c); iter++)
					{
						bool hasChild = pPS->hasChild(*iter);
						bool hasParent = pPS->hasParent(*iter);
						
						if(!hasChild && !hasParent)
						{
							// remove this sequence, it has no extensions
							pPS->removeSequence(*iter);
						}
						else if(!hasChild)
						{
							deadends.push_back(branchEnd(*iter, ANTISENSE));
						}
						else if(!hasParent)
						{
							deadends.push_back(branchEnd(*iter, SENSE));
						}
					}	
				}
			}
		}
	}	
	
	printf("seqs before trimming: %d\n", pPS->countAll());
	printf("num deadends: %d\n", deadends.size());
	
	const int MAX_DEAD_LENGTH = 2*trimNum + 5;
	
	int numBranchesRemoved = 0;
	int count = 0;
	for(std::vector<branchEnd>::iterator iter = deadends.begin(); iter != deadends.end(); iter++)
	{
		if(count % 10000 == 0)
		{
			//printf("trimming count: %d\n", count);
		}
		
		count++;
		PSequenceVector branchElements;
		
		
		extDirection dir = iter->second;
		extDirection oppositeDir = (dir == SENSE) ? ANTISENSE : SENSE;
					
		PackedSeq currSeq = iter->first;
		SeqRecord loopCheck;			
		bool stop = false;
		
		while(!stop)
		{				
			HitRecord hr = pPS->calculateExtension(currSeq, dir);
			HitRecord oppHr = pPS->calculateExtension(currSeq, oppositeDir);
				
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
				branchElements.push_back(reverseComplement(currSeq));
				//printf("stopped because of noext/ambi branch\n");
				stop = true;
			}
			else
			{
				branchElements.push_back(currSeq);
				branchElements.push_back(reverseComplement(currSeq));
					
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
				pPS->removeSequence(*bIter);
			}
		}
	}
	
	printf("seqs after trimming: %d\n", pPS->countAll());
	printf("num branches removed: %d\n", numBranchesRemoved);
}
