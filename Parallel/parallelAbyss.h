#ifndef PARALLEL_ABYSS
#define PARALLEL_ABYSS

#include "PackedSeq.h"


void doRead(int kmerSize, std::string file, int numDataNodes);
void doLoad(int sendID, int kmerSize);
void printUsage();


#endif