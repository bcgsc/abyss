#ifndef PAIREDALGORITHMS_H
#define PAIREDALGORITHMS_H

#include "Sequence.h"
#include <istream>
#include <string>
#include <vector>

struct Contig
{
	Sequence seq;
	unsigned coverage;
	Contig(Sequence seq, unsigned coverage)
		: seq(seq), coverage(coverage) { }
};

typedef std::vector<Contig> ContigVec;

typedef std::string ContigID;

namespace PairedAlgorithms
{

void readContigVec(std::string file, ContigVec& outVec);
void parseContigFromFile(std::istream& stream, ContigID& id, Sequence& seq, int& length, double& coverage);
};



#endif
