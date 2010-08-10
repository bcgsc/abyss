/** Convert various file formats to FASTQ format.
 * Written by Shaun Jackman <sjackman@bcgsc.ca>.
 * Copyright 2010 Genome Sciences Centre
 */
#include "DataLayer/Options.h"
#include "FastaReader.h"
#include "Uncompress.h"
#include <algorithm>
#include <cassert>
#include <iostream>

using namespace std;

template <class Record>
static void convert(const char* path)
{
	FastaReader in(path, FastaReader::NO_FOLD_CASE
			| FastaReader::CONVERT_QUALITY);
	for (Record record; in >> record;)
		cout << record;
	assert(in.eof());
}

int main(int argc, const char* argv[])
{
	opt::trimMasked = false;
	typedef void (*F)(const char*);
	F convertFasta = convert<FastaRecord>;
	F convertFastq = convert<FastqRecord>;
	F f = string(argv[0]).find("tofasta") != string::npos
		? convertFasta : convertFastq;
	if (argc <= 1)
		f("-");
	else
		for_each(argv + 1, argv + argc, f);
	return 0;
}
