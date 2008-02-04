#ifndef TRIMMER_H
#define TRIMMER_H

#include "Trimmer.h"
#include "Reader.h"
#include "PartitionLoader.h"
#include "PhaseSpace.h"

void ReadControlAndTrim(Config* pConfig, std::string controlFile);
void TrimByCoordinate(Config* pConfig, Coord4 trimCoord);

#endif
