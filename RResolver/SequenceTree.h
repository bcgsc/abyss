#ifndef RRESOLVER_SEQUENCETREE_H
#define RRESOLVER_SEQUENCETREE_H 1

#include "Common/Sequence.h"
#include "Contigs.h"

#include <list>
#include <vector>

std::vector<Sequence>
getTreeSequences(
    const ContigNode& start,
    const int overlap,
    const int maxLength,
    const bool forward,
    const int maxPaths);

#endif
