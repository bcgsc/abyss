#ifndef RRESOLVER_SEQUENCETREE_H
#define RRESOLVER_SEQUENCETREE_H 1

#include "Contigs.h"
#include "Common/Sequence.h"

#include <list>

std::list<Sequence>
getTreeSequences(const ContigNode& start,
                 const int overlap, const int maxLength,
                 const bool forward, const int maxPaths);

#endif
