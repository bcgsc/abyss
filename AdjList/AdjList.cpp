#include "PackedSeq.h"
#include "PairUtils.h"
#include <cerrno>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdlib.h>

using namespace std;

void generatePossibleExtensions(const PackedSeq& seq, extDirection dir, PSequenceVector& outseqs);
void readIdAndSeq(ifstream& inStream, ContigID& id, Sequence& seq);

static void assert_open(std::ifstream& f, const std::string& p)
{
	if (f.is_open())
		return;
	std::cerr << p << ": " << strerror(errno) << std::endl;
	exit(EXIT_FAILURE);
}

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
	assert_open(contigFileStream, contigFile);
	
	// Generate a k-mer -> contig lookup table for all the contig ends
	std::map<PackedSeq, ContigID> contigLUTs[2];
	
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
	  
		contigLUTs[0][seqs[0]] = id;
		contigLUTs[1][seqs[1]] = id;
	  	//std::cout << "Added S " << seqs[0].decode() << " AS " << seqs[1].decode() << std::endl;
		numAdded += 2;
		if(numAdded % 100000 == 0)
		{
			std::cout << "Added " << numAdded << std::endl;
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
					size_t lookuptable_id = 1 - idx;
					bool reverse = (compIdx == 1);
					PackedSeq testSeq;
					if(reverse)
					{
						// flip the lookup table id
						lookuptable_id = 1 - lookuptable_id;
						testSeq = reverseComplement(*iter);
					}
					else
					{
						testSeq = *iter;
					}
		  			//std::cout << id << " looking for " << testSeq.decode() << " " << lookuptable_id << std::endl;
					std::map<PackedSeq, ContigID>::iterator cLUTIter;
					cLUTIter = contigLUTs[lookuptable_id].find(testSeq);
					if(cLUTIter != contigLUTs[lookuptable_id].end())
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
