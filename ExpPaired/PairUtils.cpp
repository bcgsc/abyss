#include "PairUtils.h"
#include "Sense.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits> // for numeric_limits
#include <sstream>

using namespace std;

/** Read the distance estimates for one contig. */
istream& operator >>(istream& in, EstimateRecord& o)
{
	o.estimates[SENSE].clear();
	o.estimates[ANTISENSE].clear();

	in >> o.refID;
	in.ignore(numeric_limits<streamsize>::max(), ':');

	for (extDirection sense = SENSE; sense <= ANTISENSE; ++sense) {
		string s;
		getline(in, s, sense == SENSE ? '|' : '\n');
		istringstream ss(s);
		copy(istream_iterator<Estimate>(ss),
				istream_iterator<Estimate>(),
				back_inserter(o.estimates[sense]));
		assert(ss.eof());
	}

	return in;
}

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
		assert(serial < 1000000000);
		if (lengths.size() < serial+1)
			lengths.resize(serial+1);
		assert(lengths[serial] == 0);
		lengths[serial] = len;
	}
	assert(in.eof());
}

LinearNumKey convertContigIDToLinearNumKey(const ContigID& id)
{
	LinearNumKey key;
	key = atoi(id.c_str());
	return key;
}
