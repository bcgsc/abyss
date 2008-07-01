#include "PairUtils.h"

//
// Estimate file loaders
//

//
// Read in a single estimate from the stream
//
void readEstimateRecord(std::ifstream& stream, EstimateRecord& er)
{
	// read in the id
	ContigID cID;
	stream >> cID;
	
	er.refID = convertContigIDToNumericID(cID);
	
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

	return;
}


//
//
// Adjacency file loaders
//
//

//
//
//
void loadGraphFromAdjFile(SimpleContigGraph* pGraph, ContigLengthMap& lengthMap, std::string file)
{
	// First, load the vertices
	std::ifstream inStream(file.c_str());
	
	int numAdded = 0;
	while(!inStream.eof() && inStream.peek() != EOF)
	{
		ContigID id;
		inStream >> id;
		
		SimpleContigData data;
		data.length = lookupLength(lengthMap, id);
		NumericID numID = convertContigIDToNumericID(id);
		
		// Add the vertex to the graph
		pGraph->addVertex(numID, data);
		//pGraph->addVertex(id, empty);

		// discard the remainder of the line
		std::string discard;
		getline(inStream, discard);
		
		numAdded++;
		
		
		if(numAdded % 1000000 == 0)
		{
			printf("added %d verts\n", numAdded);
		}
	}
	
	// Now, load the edges
	
	// rewind the stream
	inStream.seekg(ios_base::beg);
	inStream.clear();
	numAdded = 0;
	while(!inStream.eof() && inStream.peek() != EOF)
	{
		ContigID id;
		inStream >> id;

		// read the adjacency info
		std::string adjRecord;
		getline(inStream, adjRecord);
		
		// parse it and add the edge
		parseAdjacencyLine(adjRecord, id, pGraph);
		
		numAdded++;
		
		if(numAdded % 1000000 == 0)
		{
			printf("added edges for %d verts\n", numAdded);
		}		
		
	}
	
	printf("end con %d %d\n", inStream.eof(), inStream.peek() != EOF);
	printf("added: %d\n", numAdded);
	size_t numVert = pGraph->getNumVertices();
	size_t numEdges = pGraph->countEdges(); // SLOW
	printf("Initial graph stats: num vert: %zu num edges: %zu\n", numVert, numEdges);	
	
}

//
//
//
void parseAdjacencyLine(std::string& adjLine, ContigID currVert, SimpleContigGraph* pGraph)
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
				ContigID adjID;
				SimpleEdgeDesc sed;
				std::stringstream recSS(record);
				recSS >> sed;
				
				// Convert the ids
				NumericID currNumID = convertContigIDToNumericID(currVert);
				NumericID adjNumID = convertContigIDToNumericID(sed.contig);
				pGraph->addEdge(currNumID, adjNumID, (extDirection)dirIdx, sed.isRC);
				//std::cout << dirIdx << " GOT ADJ: " << id << " " << reverse << "\n";
			}
		}
	}
}

//
// Length file loader
//
void loadContigLengths(std::string contigLenFile, ContigLengthMap& lengthMap)
{
	ifstream contigLenStream(contigLenFile.c_str());
	while(!contigLenStream.eof() && contigLenStream.peek() != EOF)
	{
		ContigID id;
		int len;
		contigLenStream >> id >> len;
		lengthMap[id] = len;

	}
	contigLenStream.close();
}

int lookupLength(const ContigLengthMap& lengthMap, const ContigID& id)
{
	ContigLengthMap::const_iterator iter = lengthMap.find(id);
	assert(iter != lengthMap.end());
	return iter->second;
}

//
// PDF loader
//
PDF loadPDF(std::string distCountFile, const int limit)
{
	Histogram hist;
	ifstream distFile(distCountFile.c_str());
	while(!distFile.eof() && distFile.peek() != EOF)
	{
		int value;
		int count;
		distFile >> value;
		distFile >> count;

		if(value < limit)
		{
			//std::cout << "adding " << value << " : " << count << std::endl;
			hist.addMultiplePoints(value, count);
		}
	} 

	PDF pdf(hist);
	return pdf;
}

NumericID convertContigIDToNumericID(const ContigID& id)
{
	const int nodeMult = 10000000;
	std::stringstream ss(id);
	
	int node;
	int number;
	char placeholder;
	ss >> node >> placeholder >> number;
	
	//std::cout << "id " << id << " converts to node: " << node << " num: " << number << std::endl; 
	// make sure no collisions are possible
	assert(number < nodeMult);
	return node*nodeMult + number;
}

ContigID convertNumericIDToContigID(const NumericID& id)
{
	std::stringstream ss;

	int rem = id % nodeMult;
	int pureNode = id - rem;
	int nodeID = pureNode / nodeMult;

	ss << nodeID << ":" << rem;
	ContigID ret = ss.str();
	std::cout << "ret " << ret << std::endl;
	return ret;
}
