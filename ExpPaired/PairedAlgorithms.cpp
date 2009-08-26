#include "PairedAlgorithms.h"
#include "PairUtils.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

namespace PairedAlgorithms
{
	
//
// Read contigs
//
void readContigVec(std::string file, ContigVec& outVec)
{
	std::ifstream fileHandle(file.c_str());	
	assert(fileHandle.is_open());

	while(!fileHandle.eof() && fileHandle.peek() != EOF)
	{
		ContigID strID;
		Sequence seq;
		int length;
		double coverage;
		parseContigFromFile(fileHandle, strID, seq, length, coverage);
		assert(!seq.empty());
		assert(length > 0);

		LinearNumKey numID = convertContigIDToLinearNumKey(strID);
		// The numeric ID must match the index in the vector.
		assert(numID == outVec.size());
		outVec.push_back(Contig(seq, (unsigned)coverage));
	}
	fileHandle.close();
}

void parseContigFromFile(std::istream& in,
		ContigID& id, Sequence& seq, int& length, double& coverage)
{
	assert(in.good());
	string header;
	getline(in, header);
	id.clear();
	length = 0;
	coverage = 0;
	char head;
	istringstream s(header);
   	s >> head >> id >> length >> coverage;
	assert(head == '>');

	assert(in.peek() != '>');
	seq.clear();
	while (in.peek() != '>' && !in.eof()) {
		string line;
		getline(in, line);
		seq += line;
	}
}

};
