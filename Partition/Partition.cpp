#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>

#include "Reader.h"
#include "Partition.h"
#include "IFileReader.h"
#include "FastaReader.h"
#include "PackedSeqReader.h"
#include "PackedSeqWriter.h"
#include "FastaWriter.h"
#include "PackedSeq.h"

const int MAX_FILEHANDLES = 1024;

int main(int argv, char** argc)
{
	
	if(argv < 3)
	{
		printf("usage: Partition <config file> <control file>\n");
		exit(0);
	}
	
	std::string configFile = argc[1];
	
	Config config;
	config.readConfig(configFile);	

	printf("partitioning from control\n");
	// Partition many files, as specified by the control file
	std::string controlFile = argc[2];
	PartitionFromControl(config, controlFile);
}

// Partition the input file based on the input parameters
void PartitionFile(const Config& config, std::string inputFile, int partitionSize, Coord4 start, int originalSize)
{
	printf("Partitioning %s starting at [%d %d %d %d]\n", inputFile.c_str(), start.x, start.y, start.z, start.w);
	 
	// Open a file reader
	IFileReader* pFileReader;
	
	// Check the type of file to read in
	if(inputFile.find(".fa") != -1)
	{
		//printf("opening fasta reader\n");
		pFileReader = new FastaReader(inputFile.c_str());
	}
	else if(inputFile.find(".sqb") != -1)
	{
		//printf("opening binary reader\n");
		pFileReader = new PackedSeqReader(inputFile.c_str());
	}
	else
	{
		printf("Unknown filetype: %s\n", inputFile.c_str());	
	}

	// Check that partitionSize != originalSize or this will be a waste of time
	if(partitionSize == originalSize)
	{
		printf("sub partition size matches the original size of the bins, nothing will happen, exiting (try making maxNumPartitions bigger)\n");
		exit(1);
	}	
	
	// Calculate the actual number of partitions
	int numPartitions = static_cast<int>(ceil(static_cast<double>(originalSize) / static_cast<double>(partitionSize)));
	
	printf("size: %d num: %d\n", partitionSize, numPartitions);
		
	// Sanity check, make sure we don't go over the maximum number of filehandles
	assert(pow(numPartitions,4) < MAX_FILEHANDLES);
	
	
	// Create a vector to hold the temp output file names and the true file names
	std::vector<stringPair> vecFilenames;
	
	// Create the writer objects and the directories
	
	// Place the writer pointers on a vector for linear deletion later
	// Also, create a 4D pointer space to point from a particular coordinate to a writer object
	PWV4D pWriterMap(numPartitions, PWV3D(numPartitions, PWV2D(numPartitions, PWV1D(numPartitions))));
	PWV1D pWriterVec;
	
	// Create buffers to hold the various names that will be created
	const int MAX_NAME_SIZE = 512;
	char xBuffer[MAX_NAME_SIZE];
	char yBuffer[MAX_NAME_SIZE];
	char zBuffer[MAX_NAME_SIZE];
	char wBuffer[MAX_NAME_SIZE];
	char fileBuffer[MAX_NAME_SIZE];	
	char tempFileBuffer[MAX_NAME_SIZE];			
	
	for(int x = start.x; x < start.x + originalSize; x += partitionSize)
	{	
		// Create x directory
		sprintf(xBuffer, "%s/x_%d", config.getRootDataDir().c_str(), x);
		mkdir(xBuffer, 0777);
		
		for(int y = start.y; y < start.y + originalSize; y += partitionSize)
		{
			// Create y directory
			sprintf(yBuffer, "%s/y_%d", xBuffer, y);
			mkdir(yBuffer, 0777);	
					
			for(int z = start.z; z < start.z + originalSize; z += partitionSize)
			{
				// Create z directory
				sprintf(zBuffer, "%s/z_%d", yBuffer, z);
				mkdir(zBuffer, 0777);
							
				for(int w = start.w; w < start.w + originalSize; w += partitionSize)
				{
					// Create w directory
					sprintf(wBuffer, "%s/w_%d", zBuffer, w);
					mkdir(wBuffer, 0777);

					// Create true filename
					sprintf(fileBuffer, "%s/%s", wBuffer, config.getSequenceFilename().c_str());
					
					// Create temp file to write to
					sprintf(tempFileBuffer, "%s.temp", fileBuffer);
					PackedSeqWriter* pSeqWriter = new PackedSeqWriter(tempFileBuffer, config.getSequenceLength());
					
					// add the files to the vector
					vecFilenames.push_back(stringPair(tempFileBuffer, fileBuffer));
					
					// Save the pointer for later deallocation
					pWriterVec.push_back(pSeqWriter);
					
					// Calculate the index into 4D space for this writer object
					Coord4 pos;
					pos.x = x;
					pos.y = y;
					pos.z = z;
					pos.w = w;
					
					Coord4 index = transformCoordinateToIndices(pos, start, partitionSize);
										
					assert(index.x < numPartitions && index.y < numPartitions && index.z < numPartitions && index.w < numPartitions);
					
					pWriterMap[index.x][index.y][index.z][index.w] = pSeqWriter;
				}			
			}
		}
	}
	
	int count = 0;
	//Perform the partitioning
	while(pFileReader->isGood())
	{
		// Read in a sequence
		PackedSeq pSeq = pFileReader->ReadSequence();
		
		// Ensure the sequence is reasonable
		assert(pSeq.getSequenceLength() == config.getSequenceLength());
		
		// Calculate the coordinate
		Coord4 pos = PhaseSpace::SequenceToCoord4(pSeq);
		
		// Calculate the index
		Coord4 index = transformCoordinateToIndices(pos, start, partitionSize);
		if(!(index.x < numPartitions && index.y < numPartitions && index.z < numPartitions && index.w < numPartitions))
		{
			printf("bad seq %s\n", pSeq.decode().c_str());
			printf("sequence is in incorrect bin coord: (%d %d %d %d) bin start: (%d %d %d %d)\n", pos.x, pos.y, pos.z, pos.w, start.x, start.y, start.z, start.w);
			assert(index.x < numPartitions && index.y < numPartitions && index.z < numPartitions && index.w < numPartitions);
		}
		
		pWriterMap[index.x][index.y][index.z][index.w]->WriteSequence(pSeq);
		
		if(count % 100000 == 0)
		{
			printf("partitioned %d sequences\n", count);
		}
		count++;
	}

	// Close the reader
	delete pFileReader;
	pFileReader = NULL;
		
	// Delete all the writer pointers, this will close the filehandles
	for(PWV1D::iterator iter = pWriterVec.begin(); iter != pWriterVec.end(); iter++)
	{
		delete *iter;
		*iter = 0;
	}
	
	// Move the temp files into their real location
	for(std::vector<stringPair>::iterator iter = vecFilenames.begin(); iter != vecFilenames.end(); iter++)
	{
		rename(iter->first.c_str(), iter->second.c_str());
	}
	
	// Delete the input sequence file
	//unlink(inputFile.c_str());

	return;
}

void PartitionFromControl(const Config& config, std::string controlFile)
{
	std::ifstream file(controlFile.c_str());
	
	const int MAX_LINE_LENGTH = 512;
	char buffer[MAX_LINE_LENGTH];
	
	
	while(!file.eof() && file.peek() != EOF)
	{
		file.getline(buffer, MAX_LINE_LENGTH);
		
		// Read in the values of this control line
		char seqFile[MAX_LINE_LENGTH];
		int partitionSize;
		Coord4 start;
		int originalSize;
		
		// iterate over the tokens in the line (space delimited)
		if(sscanf(buffer, "%s %d %d %d %d %d %d", seqFile, &partitionSize, &start.x, &start.y, &start.z, &start.w, &originalSize) != 7)
		{
			printf("invalid control file format\n");
			return;
		}
		
		// Perform the partition
		PartitionFile(config, seqFile, partitionSize, start, originalSize);
	}
	
	file.close();
}

Coord4 transformCoordinateToIndices(Coord4 position, Coord4 startPos, int length)
{
	Coord4 tCoord;
	tCoord.x = (position.x - startPos.x) / length;
	tCoord.y = (position.y - startPos.y) / length;
	tCoord.z = (position.z - startPos.z) / length;
	tCoord.w = (position.w - startPos.w) / length;	
	return tCoord;
} 

void test()
{
	/*
	Sequence test("ATACTAATATGCTAGCTAGTAGCTGC");
	PackedSeq encode(test);
	char shifted = encode.shiftPrepend('C');
	

	Sequence decode = encode.decode();
	printf("o: %s\n", test.c_str());
	printf("p: %s\n", decode.c_str());
	
	encode.shiftAppend(shifted);
	printf("a: %s\n", encode.decode().c_str());
	printf("shifted: %c\n", shifted);
	return 1;
	*/	
	
#if 0	
		// test as we read
		
		// check the sequence converted correctly
		assert( pSeq->decode() == newSeq);
		
		// test reverse complement
		pSeq->reverseComplement();
		assert(reverseComplement(newSeq) == pSeq->decode());
		
		// switch sequence back again and make sure it still matches
		pSeq->reverseComplement();
		
		assert(pSeq->decode() == newSeq);
		
		// test shifting
		char b = pSeq->shiftAppend('C');
		pSeq->shiftPrepend(b);
		
		assert(pSeq->decode() == newSeq);
		
		b = pSeq->shiftPrepend('T');
		pSeq->shiftAppend(b);
		
		assert(pSeq->decode() == newSeq);
		
		// Test copy constructor
		PackedSeq seq2 = *pSeq;
		
		// test operator==
		assert(seq2 == *pSeq);
#endif	
}

void createDirectoryStructure(std::string partitionDimension, int sequenceLength, int stride)
{
	// Create directories first
	/*
	char dirName[128];
	for(int i = 0; i < seqLen; i += stride)
	{
		// Directory
		sprintf(dirName, "%s/%s_%d", rootDir, partitionDimension.c_str(), i);
		mkdir(dirName, 0777);

	}
	
	printf("done directory setup\n");
	int progress = 0;
	int max = seqVector.size();
	
	for(ConstSequenceVectorIterator iter = seqVector.begin(); iter != seqVector.end(); iter++)
	{
		if(progress % 1000 == 0)
		{
			printf("progress: %d/%d\n", progress, max);
		}
		
		Coord4 c = PhaseSpace::SequenceToCoord4(*iter);
		char name[128];
		CoordToFileName(name, rootDir, c.x, c.y,c.z,c.w);
		PartitionIO file(name, FM_WRITE);
		file.WriteSequence(*iter);
		
		progress++;
	}
	*/
}

int CoordToFileName(char* buffer, const char* rootDir, int x, int y, int z, int w)
{	
	return sprintf(buffer, "%s/x_%d/y_%d/z_%d/w_%d.sqs", rootDir, x,y,z,w);
}
