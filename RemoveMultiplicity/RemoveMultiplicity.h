#ifndef REMOVE_MULTIPLICITY_H
#define REMOVE_MULTIPLICITY_H

#include <string>
#include "Config.h"
#include "CommonDefs.h"
#include "PhaseSpace.h"

void ReadControlAndRemove(Config* pConfig, std::string controlFile);
void RemoveMultiplicityByCoord(Config* pConfig, Coord4 start);

#endif
