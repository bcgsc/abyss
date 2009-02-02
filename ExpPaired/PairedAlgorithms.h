#ifndef PAIREDALGORITHMS_H
#define PAIREDALGORITHMS_H

#include "Sequence.h"
#include <fstream>
#include <string>
#include <vector>

struct Contig
{
	Sequence seq;
	bool merged;
	bool repetitive;
	bool super;
	int coverage;
};

typedef std::vector<Contig> ContigVec;

typedef std::string ContigID;

namespace PairedAlgorithms
{

void readContigVec(std::string file, ContigVec& outVec);
void parseContigFromFile(std::ifstream& stream, ContigID& id, Sequence& seq, int& length, double& coverage);
};



#endif
