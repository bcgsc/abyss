#include "PairUtils.h"
#include "Sense.h"
#include <cassert>
#include <cstdlib>
#include <fstream>
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

// Length file loader
void loadContigLengths(std::string contigLenFile, ContigLengthVec& lengthVec)
{
	ifstream contigLenStream(contigLenFile.c_str());
	assert(contigLenStream.is_open());

	while(!contigLenStream.eof() && contigLenStream.peek() != EOF)
	{
		LinearNumKey id;
		int len;		
		std::string line;
		getline(contigLenStream, line);
		
		sscanf(line.c_str(), "%d %d", &id, &len);
		
		if(id != lengthVec.size())
		{
			cerr << id << " is out of sequence (size: "
				<< lengthVec.size() << ")\n";
			assert(false);
			exit(EXIT_FAILURE);
		}

		lengthVec.push_back(len);
	}
	contigLenStream.close();
}

LinearNumKey convertContigIDToLinearNumKey(const ContigID& id)
{
	LinearNumKey key;
	key = atoi(id.c_str());
	return key;
}
