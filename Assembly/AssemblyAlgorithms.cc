#include "config.h"
#include "Common/InsOrderedMap.h"
#include <string>
#include <vector>

namespace AssemblyAlgorithms
{

/** The number of k-mer that have been eroded. */
size_t g_numEroded;

std::vector<size_t> tempCounter(16,0);
InsOrderedMap<std::string,int> tempStatMap;

void addToDb(const std::string& key, const int& value)
{
	tempStatMap.push_back(key, value);
}

};
