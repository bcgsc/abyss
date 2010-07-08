/** Convert various file formats to FASTQ format.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 * Copyright 2010 Genome Sciences Centre
 */
#include "FastaReader.h"
#include "Uncompress.h"
#include <algorithm>
#include <iostream>

using namespace std;

void convert(const char* path)
{
	FastaReader in(path, FastaReader::KEEP_N
			| FastaReader::NO_FOLD_CASE
			| FastaReader::CONVERT_QUALITY);
	for (FastqRecord fastq; in >> fastq;)
		cout << fastq;
}

int main(int argc, const char* argv[])
{
	if (argc <= 1)
		convert("-");
	else
		for_each(argv + 1, argv + argc, convert);
	return 0;
}
