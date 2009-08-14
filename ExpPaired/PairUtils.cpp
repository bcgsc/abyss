#include "PairUtils.h"
#include "Sense.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

//
// Estimate file loaders
//

//
// Read in a single estimate from the stream
//
std::istream& readEstimateRecord(std::istream& stream,
		EstimateRecord& er)
{
	er.estimates[SENSE].clear();
	er.estimates[ANTISENSE].clear();

	// read in the id
	stream >> er.refID;

	// Discard the seperator
	std::string discard;
	stream >> discard;
	
	std::string records;
	getline(stream, records);

	// convert the record to a stringstream
	std::stringstream recss(records);
	
	// Begin reading records
	size_t currIdx = 0;
	
	bool stop = false;
	while(!stop)
	{
		std::string data;
		recss >> data;
		
		if(data == "|")
		{
			currIdx = 1;
		}
		else if(data.empty())
		{
			stop = true;
		}
		else
		{
			std::stringstream dataStream(data);
			
			Estimate est;
			dataStream >> est;
			er.estimates[currIdx].push_back(est);
			//std::cout << "RECORD: " << id << " dist " << distance << " numpairs " << numPairs << std::endl;
		}
	}

	return stream;
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
