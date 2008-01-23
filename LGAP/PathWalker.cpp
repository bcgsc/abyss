#include <stdio.h>

#include <vector>
#include "Sequence.h"
#include "Reader.h"
#include "PathDriver.h"
#include "PairRecord.h"
#include "Writer.h"
#include "SeqRecord.h"
#include "Config.h"
#include "PartitionLoader.h"
#include "FastaWriter.h"
#include "FastaReader.h"

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
		
	SeqRecord extensionRecord;
	
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
												
						if(extensionRecord.contains(*iter))
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
								extensionRecord.addSequence(currSeq);				
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
					//return 0;				
				}
	printf("processed %d sequences total, %d noext %d ambiext\n", count, noext, ambiext);
	return 0;
}
