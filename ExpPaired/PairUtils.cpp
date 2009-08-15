#include "PairUtils.h"
#include "Sense.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

using namespace std;

/** Read in a single estimate from the stream. */
std::istream& readEstimateRecord(std::istream& in, EstimateRecord& o)
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
		if (serial != lengths.size()) {
			cerr << id << " is out of sequence (size: "
				<< lengths.size() << ")\n";
			exit(EXIT_FAILURE);
		}
		lengths.push_back(len);
	}
	assert(in.eof());
}

LinearNumKey convertContigIDToLinearNumKey(const ContigID& id)
{
	LinearNumKey key;
	key = atoi(id.c_str());
	return key;
}
