#ifndef PAIREDALGORITHMS_H
#define PAIREDALGORITHMS_H

#include "CommonDefs.h"
#include <fstream>
#include <string>

typedef std::vector<Contig> ContigVec;

namespace PairedAlgorithms
{

void readContigVec(std::string file, ContigVec& outVec);
void parseContigFromFile(std::ifstream& stream, ContigID& id, Sequence& seq, int& length, double& coverage);
};



#endif
