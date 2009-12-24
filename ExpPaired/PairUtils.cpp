#include "PairUtils.h"
#include "Sense.h"
#include <cassert>
#include <fstream>
#include <limits> // for numeric_limits

using namespace std;

/** Load contig lengths. */
void loadContigLengths(const string& path, ContigLengthVec& lengths)
{
	ifstream in(path.c_str());
	assert(in.is_open());

	string id;
	unsigned len;
	while (in >> id >> len) {
		in.ignore(numeric_limits<streamsize>::max(), '\n');
		LinearNumKey serial = convertContigIDToLinearNumKey(id);
		assert(serial == lengths.size());
		(void)serial;
		lengths.push_back(len);
	}
	assert(in.eof());
}
