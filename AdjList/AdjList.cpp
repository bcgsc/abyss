#include <stdio.h>
#include <math.h>
#include <iostream>
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"
#include "Options.h"
#include "FastaReader.h"
#include "PairUtils.h"

void generatePossibleExtensions(const PackedSeq& seq, extDirection dir, PSequenceVector& outseqs);
void readIdAndSeq(ifstream& inStream, ContigID& id, Sequence& seq);

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cout << "Usage: <kmer> <contig file>\n";
		return 1;
	}
	
	int kmer = atoi(argv[1]);
	std::string contigFile(argv[2]);
	
	//Load the kmer->contig map ends of the contigs
	std::cout << "File " << contigFile << " kmer " << kmer << std::endl;
	
	// Open the contig file
	ifstream contigFileStream(contigFile.c_str());
	
	// Generate a k-mer -> contig lookup table for all the contig ends
	std::map<PackedSeq, ContigID> contigLUT;
	
	int numAdded = 0;
	while(!contigFileStream.eof() && contigFileStream.peek() != EOF)
	{
		ContigID id;
		Sequence contigSequence;
		readIdAndSeq(contigFileStream, id, contigSequence);
	      
		// Generate the end->id mapping
		const unsigned numEnds = 2;
		PackedSeq seqs[numEnds];
		seqs[0] = PackedSeq(contigSequence.substr(contigSequence.length() - kmer, kmer)); //SENSE
		seqs[1] = PackedSeq(contigSequence.substr(0, kmer)); // ANTISENSE
	  
		size_t numToAdd = (seqs[0] != seqs[1]) ? 2 : 1;
	  
		numAdded += numToAdd;
		if(numAdded % 100000 == 0)
		{
			std::cout << "Added " << numAdded << std::endl;
		}
	  
	  
		for(unsigned idx = 0; idx < numToAdd; idx++)
		{	
			// insert sequences into the table
			contigLUT[seqs[idx]] = id;
		}
	} 

	// Rewind the file stream to the beginning
	contigFileStream.seekg(ios_base::beg);
	contigFileStream.clear();
	ofstream adjOutFile("AdjList.txt");

	int numVerts = 0;
	int numEdges = 0;

	// Build the edges and write them out
	while(!contigFileStream.eof() && contigFileStream.peek() != EOF)
	{
		ContigID id;
		Sequence contigSequence;
		readIdAndSeq(contigFileStream, id, contigSequence);
		adjOutFile << id << " ";
		// Generate edges to/from this node
  
		// Since two contigs are not necessarily built from the same strand, two contigs can both have OUT nodes pointing to each other
		// this situation will get cleaned up when the links are resolved/merged
		const unsigned numEnds = 2;
		PackedSeq seqs[numEnds];
		seqs[0] = PackedSeq(contigSequence.substr(contigSequence.length() - kmer, kmer)); //SENSE
		seqs[1] = PackedSeq(contigSequence.substr(0, kmer)); // ANTISENSE
		  
		ExtensionRecord extRec;
		  
		for(unsigned idx = 0; idx < numEnds; idx++)
		{
			std::vector<SimpleEdgeDesc> edges;
			PackedSeq& currSeq = seqs[idx];
			extDirection dir;
			dir = (idx == 0) ? SENSE : ANTISENSE;
			  
			// Generate the links
			PSequenceVector extensions;
			generatePossibleExtensions(currSeq, dir, extensions);
		  
			for(PSequenceVector::iterator iter = extensions.begin(); iter != extensions.end(); ++iter)
			{
				// Get the contig this sequence maps to
				for(size_t compIdx = 0; compIdx <= 1; ++compIdx)
				{
					bool reverse = (compIdx == 1);
					PackedSeq testSeq;
					if(reverse)
					{
						testSeq = reverseComplement(*iter);
					}
					else
					{
						testSeq = *iter;
					}
		  
					std::map<PackedSeq, ContigID>::iterator cLUTIter;
					cLUTIter = contigLUT.find(testSeq);
					if(cLUTIter != contigLUT.end())
					{						
						// Store the edge
						SimpleEdgeDesc ed;
						ed.contig = cLUTIter->second;
						ed.isRC = reverse;
						edges.push_back(ed);
					}
				}
			}
			// Print the edges
			adjOutFile << "[ ";
			std::copy(edges.begin(), edges.end(), std::ostream_iterator<SimpleEdgeDesc>(adjOutFile, " "));
			adjOutFile << "] ";
			numEdges += edges.size();
		}
		adjOutFile << "\n";
		numVerts++;
	}

	printf("Found %d edges for %d verts\n", numEdges, numVerts);
  
	contigFileStream.close();
	adjOutFile.close();
} 

void readIdAndSeq(ifstream& inStream, ContigID& id, Sequence& seq)
{

  // Read the contig id and the sequence
  std::string tempName;
  inStream >> tempName;
      
  // Chop the first character to produce the name
  id = tempName.substr(1);
  
  // Read to the end of the line
  std::string discard;
  getline(inStream, discard);
  
  // Read in the sequence
  getline(inStream, seq);
  
  //std::cout << "ID: " << id << " seq: " << seq << std::endl;
}

void generatePossibleExtensions(const PackedSeq& seq, extDirection dir, PSequenceVector& outseqs)
{
  PackedSeq extSeq(seq);
  extSeq.rotate(dir, 'A');
  
  // Check for the existance of the 4 possible extensions
  for(int i  = 0; i < NUM_BASES; i++)
    {
      char currBase = BASES[i];
      extSeq.setLastBase(dir, currBase);
      outseqs.push_back(extSeq);
    }
}
