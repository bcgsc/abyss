#ifndef ASSEMBLY_BRANCHRECORD_H
#define ASSEMBLY_BRANCHRECORD_H 1

/** Generate the sequence of this contig. */
template <typename It, typename OutIt>
void branchRecordToStr(It it, It last, OutIt out)
{
	assert(it < last);
	std::string k0 = it->first.str();
	std::copy(k0.begin(), k0.end(), out);
	out += k0.length();
	++it;
	for (; it != last; ++it) {
		*out = it->first.getLastBaseChar();
		++out;
	}
}

#include "Assembly/BranchRecordBase.h"

#endif
