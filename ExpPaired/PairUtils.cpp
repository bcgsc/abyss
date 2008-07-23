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
void loadGraphFromAdjFile(SimpleContigGraph* pGraph,  std::string& lengthFile, std::string adjFile)
{
	
	// Load the lengths temporarily
	ContigLengthVec* pLengthVec = new ContigLengthVec();
	loadContigLengths(lengthFile, *pLengthVec);
	
	// First, load the vertices
	std::ifstream inStream(adjFile.c_str());
	
	int numAdded = 0;
	while(!inStream.eof() && inStream.peek() != EOF)
	{
		LinearNumKey id;
		inStream >> id;
		
		SimpleContigData data;
		data.length = lookupLength(*pLengthVec, id);
		
		// Add the vertex to the graph
		pGraph->addVertex(id, data);
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
	
	// Delete the lengths to free up space
	delete pLengthVec;
	pLengthVec = NULL;
	
	// Now, load the edges
	
	// rewind the stream
	inStream.seekg(ios_base::beg);
	inStream.clear();
	numAdded = 0;
	while(!inStream.eof() && inStream.peek() != EOF)
	{
		LinearNumKey id;
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
	std::cout << "Loading lengths\n";
	ifstream contigLenStream(contigLenFile.c_str());
	while(!contigLenStream.eof() && contigLenStream.peek() != EOF)
	{
		LinearNumKey id;
		int len;		
		std::string line;
		getline(contigLenStream, line);
		
		sscanf(line.c_str(), "%d %d", &id, &len);
		
		if(id != lengthVec.size())
		{
			std::cout << id << " is out of sequence (size: " << lengthVec.size() << ")\n";
			assert(false);
		}

		lengthVec.push_back(len);

	}
	contigLenStream.close();
	std::cout << "Done loading lengths\n";
}

int lookupLength(const ContigLengthVec& lengthVec, const LinearNumKey& id)
{
	assert(id < lengthVec.size());
	return lengthVec[id];
}

//
// PDF loader
//
PDF loadPDF(std::string distCountFile)
{
	Histogram hist;
	ifstream distFile(distCountFile.c_str());
	while(!distFile.eof() && distFile.peek() != EOF)
	{
		int value;
		int count;
		distFile >> value;
		distFile >> count;
		hist.addMultiplePoints(value, count);
	} 

	PDF pdf(hist);
	return pdf;
}

LinearNumKey convertContigIDToLinearNumKey(const ContigID& id)
{
	LinearNumKey key;
	key = atoi(id.c_str());
	return key;
}
