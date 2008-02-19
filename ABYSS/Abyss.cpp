#include <stdio.h>

#include <vector>
#include <stdio.h>
#include <deque>
#include <iostream>
#include <fstream>
#include "Abyss.h"
#include "CommonUtils.h"
#include "SimpleSequenceSpace.h"

ofstream branchLog("branchLog.txt");

int main(int argc, char** argv)
{	

	if(argc < 5 || argv[1] == "--help")
	{
		printUsage();
		exit(1);
	}

	std::string fastaFile = argv[1];
	int readLen = atoi(argv[2]);
	int kmerSize = atoi(argv[3]);
	int numTrims = atoi(argv[4]);

	// Set up coordinate space
	int maxC = kmerSize - 1;
	Coord4 minCoords = {0, 0, 0, 0};
	Coord4 maxCoords = {maxC, maxC, maxC, maxC};
	
	// Load the phase space
	SimpleSequenceSpace* pSS = new SimpleSequenceSpace(kmerSize, minCoords, maxCoords);
	
	// Open file for read
	FastaReader* reader = new FastaReader(fastaFile.c_str());
	
	// Load phase space
	int count = 0;
	
	PSequenceVector seqs;
	
	while(reader->isGood())
	{
		PackedSeq seq = reader->ReadSequence();		
		assert(kmerSize <= seq.getSequenceLength());
		
		for(int i = 0; i < seq.getSequenceLength() - kmerSize  + 1; i++)
		{
			PackedSeq sub = seq.subseq(i, kmerSize);
			pSS->addSequence(sub, true);

			count++;
			
			if(count % 1000000 == 0)
			{
				printf("loaded %d sequences\n", count);
			}
			
			seqs.push_back(sub);
		}
	}

	printf("total sequences: %d\n", pSS->countAll());
	
	printf("finalizing\n");
	pSS->finalizeBins();	

	printf("generating adjacency info\n");
	pSS->generateAdjacency();
	
	count = 0;

	//for(SequenceCollectionIter iter = pSS->getStartIter(); iter != pSS->getEndIter(); iter++)
	//{
		//assert(pSS->checkForSequence(*iter));
	//}

	
	// close the reader	
	delete reader;

	// trim the reads
	outputBranchSizes(pSS, minCoords, maxCoords);
	//for(int i = 0; i < numTrims; i++)
	int start = 2;
	int step = 2;
	int maxBranch =  4 * (readLen - kmerSize + 1);
	
	
	while(start <= maxBranch)
	{
		//trimSequences(pSS, minCoords, maxCoords);
		trimSequences2(pSS, minCoords, maxCoords, start);
		//trimSequences3(pSS, minCoords, maxCoords, i+1, 3 * (readLen - kmerSize + 1));
		start << 1;
	}
	
	// Now trim at the max branch length
	for(int i = 0; i < 2; i++)
	{
		trimSequences2(pSS, minCoords, maxCoords, maxBranch);
	}
	
	// finally, trim at the min contig length
	trimSequences2(pSS, minCoords, maxCoords, 100);
	
	outputBranchSizes(pSS, minCoords, maxCoords);
	
	printf("outputting trimmed reads\n");
	outputSequences("trimmed.fa", pSS, minCoords, maxCoords);
	
	printf("assembling sequences\n");
	assemble(pSS, minCoords, maxCoords);

	delete pSS;
	return 0;
}

void outputSequences(const char* filename, SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord)
{
	FastaWriter writer(filename);

					
	for(SequenceCollectionIter iter = pSS->getStartIter(); iter != pSS->getEndIter(); iter++)
	{
		if(!pSS->checkSequenceFlag(*iter, SF_DELETE))
		{
			writer.WriteSequence(*iter);
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

void assemble(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord)
{	
	// create file writer
	FastaWriter writer("contigs.fa");
	printf("starting assembly\n");
	for(SequenceCollectionIter iter = pSS->getStartIter(); iter != pSS->getEndIter(); iter++)
	{
		if(!pSS->checkSequenceFlag(*iter, SF_SEEN) && !pSS->checkSequenceFlag(*iter, SF_DELETE))
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
					pSS->markSequence(currSeq, SF_SEEN);
									
					HitRecord hr = pSS->calculateExtension(currSeq, dir);
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
				writer.WriteSequence(contig);
			}
		}		
	}
	
	printf("noext: %d, ambi: %d\n", noext, ambiext);
}

void assemble2(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord)
{	
	std::list<branchEnd> seqStarts;
			
	for(SequenceCollectionIter iter = pSS->getStartIter(); iter != pSS->getEndIter(); iter++)
	{
		bool hasChild = pSS->hasChild(*iter);
		bool hasParent = pSS->hasParent(*iter);
						
		if(!hasChild)
		{
			seqStarts.push_back(branchEnd(*iter, ANTISENSE));
		}
		else if(!hasParent)
		{
			seqStarts.push_back(branchEnd(*iter, SENSE));
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
			
		if(!pSS->checkSequenceFlag(currSeq, SF_SEEN))
		{
			printf("starting from %s\n", currSeq.decode().c_str());
			
			// the record of extensions			
			PSequenceVector extensions;

			bool stop = false;
			SeqRecord loopCheck;			
				
			while(!stop)
			{
					
				// Mark the flag for the selected sequence
				pSS->markSequence(currSeq, SF_SEEN);
				pSS->markSequence(reverseComplement(currSeq), SF_SEEN);
									
				HitRecord hr = pSS->calculateExtension(currSeq, dir);
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

Sequence assembleSequence(SimpleSequenceSpace* pSS, SequenceCollectionIter sequenceIter)
{		
	Sequence contig;
	return contig;
}

void trimSequences(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord)
{				
	for(SequenceCollectionIter iter = pSS->getStartIter(); iter != pSS->getEndIter(); iter++)
	{
		if(!iter->isFlagSet(SF_DELETE))
		{
			if(!(pSS->hasChild(*iter) && pSS->hasParent(*iter)))
			{
				pSS->removeSequence(*iter);
			}
		}
	}	
}

void trimSequences2(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord, int maxBranchCull)
{
	const int MAX_DEAD_LENGTH = maxBranchCull;
	printf("trimming max branch: %d\n", maxBranchCull);	
	int numBranchesRemoved = 0;
	int count = 0;
	for(SequenceCollectionIter iter = pSS->getStartIter(); iter != pSS->getEndIter(); iter++)
	{
		if(iter->isFlagSet(SF_DELETE))
		{
			continue;
		}
				
		bool hasChild = pSS->hasChild(*iter);
		bool hasParent = pSS->hasParent(*iter);
		
		extDirection dir;
		if(!hasChild && !hasParent)
		{
			// remove this sequence, it has no extensions
			pSS->removeSequence(*iter);
			continue;
		}
		else if(!hasChild)
		{
			dir = ANTISENSE;
		}
		else if(!hasParent)
		{
			dir = SENSE;
		}
		else
		{
			continue;	
		}

		count++;
		PSequenceVector branchElements;
		
		extDirection oppositeDir = oppositeDirection(dir);
					
		PackedSeq currSeq = *iter;
		SeqRecord loopCheck;			
		bool stop = false;
		
		while(!stop)
		{				
			HitRecord hr = pSS->calculateExtension(currSeq, dir);
			HitRecord oppHr = pSS->calculateExtension(currSeq, oppositeDir);
				
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
				pSS->removeSequence(*bIter);
			}
		}
	}
	
	printf("seqs after trimming: %d\n", pSS->countAll());
	printf("num branches removed: %d\n", numBranchesRemoved);
}


void outputBranchSizes(SimpleSequenceSpace* pSS, Coord4 minCoord, Coord4 maxCoord)
{	 
	static int trimNum = 0;
	//printf("seqs before trimming: %d\n", pSS->countAll());
	int numBranchesRemoved = 0;
	int count = 0;
	for(SequenceCollectionIter iter = pSS->getStartIter(); iter != pSS->getEndIter(); iter++)
	{	
		if(iter->isFlagSet(SF_DELETE))
		{
			continue;
		}
				
		// does this sequence have an parent/child extension?
		bool hasChild = pSS->hasChild(*iter);
		bool hasParent = pSS->hasParent(*iter);
		extDirection dir;
		
		if(!hasChild && !hasParent)
		{
			// remove this sequence, it has no extensions
			pSS->removeSequence(*iter);
			continue;
		}
		else if(!hasChild)
		{
			dir = ANTISENSE;
		}
		else if(!hasParent)
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
			HitRecord hr = pSS->calculateExtension(currSeq, dir);
			HitRecord oppHr = pSS->calculateExtension(currSeq, oppositeDir);
				
			if(oppHr.getNumHits() > 1)
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
		branchLog << trimNum << " " << branchElements.size() << endl;
	}
	
	trimNum++;
}



void printUsage()
{
	printf("usage: ABYSS <reads fasta file> <max read length> <kmer size> <number of trimming steps>\n");	
}
