#ifndef PARTITION_H
#define PARTITION_H

#include <fstream>
#include "IFileWriter.h"
#include "PackedSeqWriter.h"
#include "Config.h"

typedef std::pair<std::string, std::string> stringPair;

typedef std::vector<IFileWriter*> PWV1D;
typedef std::vector<PWV1D> PWV2D;
typedef std::vector<PWV2D> PWV3D;
typedef std::vector<PWV3D> PWV4D;

void PartitionFile(const Config& config, std::string inputFile, int partitionSize, Coord4 start, int originalSize);
void PartitionFromControl(const Config& config, std::string controlFile);

int CoordToFileName(char* buffer, const char* rootDir, int x, int y, int z, int w);
Coord4 transformCoordinateToIndices(Coord4 position, Coord4 startPos, int length);
void createDirectoryStructure(std::string partitionDimension, int sequenceLength, int stride);

#endif

