#include "AssemblyAlgorithms.h"
#include "Assembly/Options.h"
#include "Common/Options.h"
#include "FastaReader.h"
#include "FastaWriter.h"
#include "Histogram.h"
#include "IOUtil.h"
#include "Log.h"
#include "SequenceCollection.h"
#include "StringUtil.h"
#include "Timer.h"
#include <algorithm>
#include <cctype>
#include <climits> // for UINT_MAX
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace AssemblyAlgorithms
{

/** The number of k-mer that have been eroded. */
size_t g_numEroded;

#if _SQL
std::vector<size_t> tempCounter(16,0);
InsOrderedMap<std::string,int> tempStatMap;

void addToDb(const std::string& key, const int& value)
{
	tempStatMap.push_back(key, value);
}
#endif

};
