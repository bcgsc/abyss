#include <stdio.h>

#include <vector>
#include <stdio.h>
#include "PathWalker.h"


#if 0
int main(int argv, char** argc)
{
	if(argv < 2)
	{
		printf("usage: Trimmer <config file>\n");
		exit(1);
	}
	
	std::string configFile = argc[1];
	
	// Read in the config
	Config config;
	config.readConfig(configFile);	
	
	// Calculate the coordinate extents, all the sequences in this range will be loaded. these coordinates are INCLUSIVE
	Coord4 minCoords;
	Coord4 maxCoords;
	Coord4 trimCoord = {2, 8, 10, 16};
	
	// get the step size
	int stepSize = config.getUnitSize();
	
	// Compute the coordinate extents
	minCoords.x = max(trimCoord.x - 2, 0);
	maxCoords.x = min(trimCoord.x + stepSize + 2, config.getSequenceLength() -1);
	
	minCoords.y = max(trimCoord.y - 1, 0);
	maxCoords.y = min(trimCoord.y + stepSize, config.getSequenceLength() - 1);
	
	minCoords.z = max(trimCoord.z - 1, 0);
	maxCoords.z = min(trimCoord.z + stepSize, config.getSequenceLength() - 1);
	
	minCoords.w = max(trimCoord.w - 1, 0);
	maxCoords.w = min(trimCoord.w + stepSize, config.getSequenceLength() - 1);	

	printf("trim coord: (%d, %d, %d, %d)\n", trimCoord.x, trimCoord.y, trimCoord.z, trimCoord.w);
	printf("min coord: (%d, %d, %d, %d)\n", minCoords.x, minCoords.y, minCoords.z, minCoords.w);
	printf("max coord: (%d, %d, %d, %d)\n", maxCoords.x, maxCoords.y, maxCoords.z, maxCoords.w);
	
	PartitionLoader pl(&config);
	PhaseSpace* pPS = pl.CreateAndLoadPhaseSpace(minCoords, maxCoords);	
	
	// Create the output file
	std::string outFile = PartitionLoader::Coord4ToPartitionFile(&config, trimCoord);
	std::string tempFile = outFile + config.getTempFileExtension();	
	
	PackedSeqWriter writer(tempFile.c_str(), config.getSequenceLength());
	printf("opened %s for write\n", tempFile.c_str());
	
	Coord4 size;
	size.x = min(stepSize, config.getSequenceLength() - trimCoord.x);
	size.y = min(stepSize, config.getSequenceLength() - trimCoord.y);
	size.z = min(stepSize, config.getSequenceLength() - trimCoord.z);
	size.w = min(stepSize, config.getSequenceLength() - trimCoord.w);
	 
	int noext = 0;
	int ambiext = 0;
	int count = 0;
	int boundary = 0;
		 
	for(int x = 0; x < size.x; x++)
		for(int y = 0; y < size.y; y++)
			for(int z = 0; z < size.z; z++)
				for(int w = 0; w < size.w; w++)
				{
					Coord4 c;
					c.x = trimCoord.x + x;
					c.y = trimCoord.y + y;
					c.z = trimCoord.z + z;
					c.w = trimCoord.w + w;
					
					for(PhaseSpaceBinIter iter = pPS->getStartIter(c); iter != pPS->getEndIter(c); iter++)
					{											
						if(iter->isFlagSet(SF_SEEN))
						{
							continue;
						}
						
						// To hold the extension path
						PSequenceVector extensions[2];
						
						for(int i = 0; i <= 1; i++)
						{
							bool stop = false;
							
							extDirection dir = (i == 0) ? SENSE : ANTISENSE;
							PackedSeq currSeq = *iter;
							SeqRecord loopCheck;
							
							while(!stop)
							{
								HitRecord hr = pPS->calculateExtension(currSeq, dir);
								if(hr.getNumHits() == 0)
								{
									// no ext
									stop = true;
									noext++;
								}
								else if(hr.getNumHits() == 1)
								{
									// move curr sequence 
									currSeq = hr.getFirstHit();				
				
									extensions[i].push_back(currSeq);							
									
									Coord4 currCoord = PhaseSpace::SequenceToCoord4(currSeq);	
									// check if we should stop because of loops or hit a boundary
									
									if(loopCheck.contains(currSeq) || !isCoordInternal(currCoord, trimCoord, size))
									{
										stop = true;
									}
									else
									{
										pPS->markSequence(currSeq, SF_DELETE);
										pPS->markSequence(reverseComplement(currSeq), SF_DELETE);
									}
								}
								else
								{
									// ambi ext
									stop = true;
									ambiext++;
								}
								
								loopCheck.addSequence(currSeq);
								
								// Mark the flag for the selected sequence
								pPS->markSequence(currSeq, SF_SEEN);
								pPS->markSequence(reverseComplement(currSeq), SF_SEEN);
							}
						}
						
						Sequence contig;
						contig.reserve(iter->getSequenceLength() + extensions[0].size() + extensions[1].size());
						
						// output the contig
						
						// output all the antisense extensions
						for(PSequenceVector::reverse_iterator asIter = extensions[1].rbegin(); asIter != extensions[1].rend(); asIter++)
						{
							contig.append(1, asIter->getFirstBase());
							
						}
						
						// output the current sequence itself
						contig.append(iter->decode());
						
						// output the sense extensions
						for(PSequenceVector::iterator sIter = extensions[0].begin(); sIter != extensions[0].end(); sIter++)
						{
							contig.append(1, sIter->getLastBase());
						}
						
						writer.WriteSequence(contig);
					}					
				}
	
	printf("processed %d sequences total, %d noext %d ambiext boundary %d\n", count, noext, ambiext, boundary);
	
	int cull = 0;
	int notcull = 0;
	for(int x = 0; x < size.x; x++)
		for(int y = 0; y < size.y; y++)
			for(int z = 0; z < size.z; z++)
				for(int w = 0; w < size.w; w++)
				{
					Coord4 c;
					c.x = trimCoord.x + x;
					c.y = trimCoord.y + y;
					c.z = trimCoord.z + z;
					c.w = trimCoord.w + w;
					
					for(PhaseSpaceBinIter iter = pPS->getStartIter(c); iter != pPS->getEndIter(c); iter++)
					{	
						if(iter->isFlagSet(SF_DELETE))
						{
							cull++;
						}
						else
						{
							notcull++;	
						}
					}
				}
				
	printf("%d cull %d not cull\n", cull, notcull);
	delete pPS;
	pPS = 0;	
}
#endif

bool isCoordInternal(Coord4 c, Coord4 start, Coord4 size)
{
	if(c.x >= start.x && c.x < start.x + size.x &&
	   c.y >= start.y && c.y < start.y + size.y &&
	   c.z >= start.z && c.z < start.z + size.z &&
	   c.w >= start.w && c.w < start.w + size.w)
	{
		return true;
	}
	else
	{
		return false;
	}
	
}

int main(int argv, char** argc)
{	
	if(argv < 2)
	{
		printf("usage: LGAP <fasta file> <config file>\n");
		exit(1);
	}
	
	std::string fastaFile = argc[1];
	std::string configFile = argc[2];
	
	Config config;
	config.readConfig(configFile);
	
	Coord4 minCoords;
	Coord4 maxCoords;
	
	int maxC = config.getSequenceLength() - 1;
	// Compute the coordinate extents
	minCoords.x = 0;
	maxCoords.x = maxC;
	
	minCoords.y = 0;
	maxCoords.y = maxC;
	
	minCoords.z = 0;
	maxCoords.z = maxC;
	
	minCoords.w = 0;
	maxCoords.w = maxC;	
	
	
	// Load the phase space
	PhaseSpace* pPS = new PhaseSpace(config.getSequenceLength(), minCoords, maxCoords);
	
	// Open file for read
	FastaReader* reader = new FastaReader(fastaFile.c_str());
	
	// Load phase space
	int count = 0;
	
	while(reader->isGood())
	{
		PackedSeq seq = reader->ReadSequence();
		pPS->addSequence(seq, true);
		count++;
	}
	
	delete reader;
	
	printf("done phase space load (%d sequences)\n", count);
	
	pPS->finalizeBins(minCoords, maxCoords);
	
	// create file writer
	FastaWriter writer("contigs.fa");
	
	printf("done load\n");

	int noext = 0;
	int ambiext = 0;
	
	count = 0;
	for(int x = 0; x < maxC; x++)
		for(int y = 0; y < maxC; y++)
			for(int z = 0; z < maxC; z++)
				for(int w = 0; w < maxC; w++)
				{
					Coord4 c;
					c.x = x;
					c.y = y;
					c.z = z;
					c.w = w;
					
					for(PhaseSpaceBinIter iter = pPS->getStartIter(c); iter != pPS->getEndIter(c); iter++)
					{
						count++;
						if(count % 100000 == 0)
						{
							printf("processed %d sequences, %d noext %d ambiext\n", count, noext, ambiext);
						}
											
												
						if(iter->isFlagSet(SF_SEEN))
						{
							continue;
						}
						
						PSequenceVector extensions[2];
							
						int stage = 0;
						for(int i = 0; i <= 1; i++)
						{
							bool stop = false;
							
							extDirection dir = (i == 0) ? SENSE : ANTISENSE;
							PackedSeq currSeq = *iter;
							SeqRecord loopCheck;
							
							while(!stop)
							{
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
									currSeq = hr.getFirstHit();
									
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
								
								// Mark the flag for the selected sequence
								pPS->markSequence(currSeq, SF_SEEN);
								pPS->markSequence(reverseComplement(currSeq), SF_SEEN);
							}
						}
						
						Sequence contig;
						contig.reserve(iter->getSequenceLength() + extensions[0].size() + extensions[1].size());
						
						// output the contig
						
						// output all the antisense extensions
						for(PSequenceVector::reverse_iterator asIter = extensions[1].rbegin(); asIter != extensions[1].rend(); asIter++)
						{
							contig.append(1, asIter->getFirstBase());
							
						}
						
						// output the current sequence itself
						contig.append(iter->decode());
						
						// output the sense extensions
						for(PSequenceVector::iterator sIter = extensions[0].begin(); sIter != extensions[0].end(); sIter++)
						{
							contig.append(1, sIter->getLastBase());
						}
						
						if(contig.length() >= 100)
						{
							writer.WriteSequence(contig);
						}
					}
					//return 0;				
				}
	printf("processed %d sequences total, %d noext %d ambiext\n", count, noext, ambiext);
	return 0;
}
