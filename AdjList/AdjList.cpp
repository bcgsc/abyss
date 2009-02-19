#include "PackedSeq.h"
#include "PairUtils.h"
#include <cassert>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

static int opt_verbose = 0;

void generatePossibleExtensions(const PackedSeq& seq, extDirection dir, PSequenceVector& outseqs);
void readIdAndSeq(ifstream& inStream, ContigID& id, Sequence& seq);

static void assert_open(std::ifstream& f, const std::string& p)
{
	if (f.is_open())
		return;
	std::cerr << p << ": " << strerror(errno) << std::endl;
	exit(EXIT_FAILURE);
}

/** A contig ID and its end sequences. */
struct ContigEndSeq {
	ContigID id;
	PackedSeq l;
	PackedSeq r;
	ContigEndSeq(const ContigID& id,
			const PackedSeq& l, const PackedSeq& r)
		: id(id), l(l), r(r) { }
};

int main(int argc, char** argv)
{
	if(argc < 3)
	{
		std::cerr << "Usage: <kmer> <contig file>\n";
		return 1;
	}

	int kmer = atoi(argv[1]);
	string contigPath(argv[2]);

	ifstream contigFileStream(contigPath.c_str());
	assert_open(contigFileStream, contigPath);

	// Generate a k-mer -> contig lookup table for all the contig ends
	vector<ContigEndSeq> contigs;
	std::map<PackedSeq, ContigID> contigLUTs[2];

	while(!contigFileStream.eof() && contigFileStream.peek() != EOF)
	{
		ContigID id;
		Sequence seq;
		readIdAndSeq(contigFileStream, id, seq);
		assert(contigFileStream.good());
		PackedSeq seql = seq.substr(seq.length() - kmer, kmer);
		PackedSeq seqr = seq.substr(0, kmer);
		contigs.push_back(ContigEndSeq(id, seql, seqr));
		contigLUTs[0][seql] = id;
		contigLUTs[1][seqr] = id;
	}
	assert(contigFileStream.eof());
	contigFileStream.close();

	ostream& out = cout;
	int numVerts = 0;
	int numEdges = 0;
	for (vector<ContigEndSeq>::const_iterator i = contigs.begin();
			i != contigs.end(); ++i) {
		const ContigID& id = i->id;
		const PackedSeq seqs[2] = { i->l, i->r };

		out << id;

		const unsigned numEnds = 2;
		for(unsigned idx = 0; idx < numEnds; idx++)
		{
			std::vector<SimpleEdgeDesc> edges;
			const PackedSeq& currSeq = seqs[idx];
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
			out << " [ ";
			std::copy(edges.begin(), edges.end(), std::ostream_iterator<SimpleEdgeDesc>(out, " "));
			out << ']';
			numEdges += edges.size();
		}
		out << '\n';
		numVerts++;
	}

	if (opt_verbose > 0)
		cerr << "vertices: " << numVerts << " "
			"edges: " << numEdges << endl;
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
