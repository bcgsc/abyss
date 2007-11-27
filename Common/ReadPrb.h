#ifndef READPRB_H
#define READPRB_H

#include <vector>
#include "Prb.h"

typedef std::vector<Prb> ReadPrb;

double calculatePError(const ReadPrb& readPrb);

#endif
