#include "PairUtils.h"
#include <cassert>
#include <fstream>

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


//
//
// Adjacency file loaders
//
//

//
//
//
void loadGraphFromAdjFile(SimpleContigGraph* pGraph,  std::string& lengthFile, std::string adjFile)
{
	// Load the lengths temporarily
	ContigLengthVec* pLengthVec = new ContigLengthVec();
	loadContigLengths(lengthFile, *pLengthVec);

	// First, load the vertices
	std::ifstream inStream(adjFile.c_str());
	assert(inStream.is_open());

	int numAdded = 0;
	LinearNumKey id;
	std::string adjRecord;
	while (inStream >> id
			&& getline(inStream, adjRecord)) {
		SimpleContigData data;
		data.length = pLengthVec->at(id);
		pGraph->addVertex(id, data);

		numAdded++;
		if (numAdded % 1000000 == 0)
			printf("added %d verts\n", numAdded);
	}
	assert(inStream.eof());

	// Delete the lengths to free up space
	delete pLengthVec;
	pLengthVec = NULL;

	// Now, load the edges
	inStream.clear();
	inStream.seekg(ios_base::beg);
	numAdded = 0;
	while (inStream >> id
			&& getline(inStream, adjRecord)) {
		parseAdjacencyLine(adjRecord, id, pGraph);

		numAdded++;
		if (numAdded % 1000000 == 0)
			printf("added edges for %d verts\n", numAdded);
	}
	assert(inStream.eof());

	size_t numVert = pGraph->getNumVertices();
	size_t numEdges = pGraph->countEdges(); // SLOW
	printf("Initial graph stats: num vert: %zu num edges: %zu\n",
			numVert, numEdges);
}

//
//
//
void parseAdjacencyLine(std::string& adjLine, LinearNumKey currVert, SimpleContigGraph* pGraph)
{
	//std::cout << "ADJ RECORD: " << adjLine << std::endl;
	
	// convert to string stream
	std::stringstream ss(adjLine);
	for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
	{
		// Extract the opening bracket
		std::string bracket;
		ss >> bracket;
		
		// Begin extracting adjacency records
		bool done = false;
		while(!done)
		{
			// extract the record
			std::string record;
			ss >> record;
			
			// check if the record is valid or we've hit the end
			if(record == "]")
			{
				done = true;
			}
			else
			{
				LinearNumKey adjID;
				SimpleEdgeDesc sed;
				std::stringstream recSS(record);
				recSS >> sed;
				
				adjID = convertContigIDToLinearNumKey(sed.contig);
				
				// Convert the ids
				pGraph->addEdge(currVert, adjID, (extDirection)dirIdx, sed.isRC);
				//std::cout << dirIdx << " GOT ADJ: " << id << " " << reverse << "\n";
			}
		}
	}
}

//
// Length file loader
//
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

//
// Hist loader
//
Histogram loadHist(std::string distCountFile)
{
	ifstream distFile(distCountFile.c_str());
	assert(distFile.is_open());

	Histogram hist;
	int value;
	unsigned count;
	while (distFile >> value >> count)
		hist.addMultiplePoints(value, count);
	assert(distFile.eof());
	return hist;
}

LinearNumKey convertContigIDToLinearNumKey(const ContigID& id)
{
	LinearNumKey key;
	key = atoi(id.c_str());
	return key;
}
