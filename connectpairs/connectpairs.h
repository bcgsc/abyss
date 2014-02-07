#ifndef CONNECTPAIRS_H
#define CONNECTPAIRS_H
#include "DataLayer/FastaInterleave.h"

/** Uppercase only bases that are present in original reads.
 *  @return number of mis-matching bases. */
static inline unsigned maskNew(const FastqRecord& read1, const FastqRecord& read2,
		FastaRecord& merged, int mask = 0)
{
	Sequence r1 = read1.seq, r2 = reverseComplement(read2.seq);
	if (mask) {
		transform(r1.begin(), r1.end(), r1.begin(), ::tolower);
		transform(r2.begin(), r2.end(), r2.begin(), ::tolower);
		transform(merged.seq.begin(), merged.seq.end(), merged.seq.begin(),
				::tolower);
	}
	unsigned mismatches = 0;
	for (unsigned i = 0; i < r1.size(); i++) {
		assert(i < merged.seq.size());
		if (r1[i] == merged.seq[i])
			merged.seq[i] = toupper(r1[i]);
		else
			mismatches++;
	}
	for (unsigned i = 0; i < r2.size(); i++) {
		assert(r2.size() <= merged.seq.size());
		unsigned merged_loc = i + merged.seq.size() - r2.size();
		if (r2[i] == merged.seq[merged_loc])
			merged.seq[merged_loc] = toupper(r2[i]);
		else
			mismatches++;
	}
	return mismatches;
}

#endif
