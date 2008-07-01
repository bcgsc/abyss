#ifndef LOADERS_H
#define LOADERS_H

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include "SequenceCollectionHash.h"
#include "AssemblyAlgorithms.h"
#include "Options.h"
#include "FastaReader.h"
#include "Stats.h"
#include "DirectedGraph.h"

const int nodeMult = 10000000;

struct SimpleContigData
{
	int length;
};

typedef std::map<ContigID, int> ContigLengthMap;
typedef int NumericID;
typedef DirectedGraph<NumericID, SimpleContigData> SimpleContigGraph;

// FUNCTIONS
NumericID convertContigIDToNumericID(const ContigID& id);
ContigID convertNumericIDToContigID(const NumericID& id);

// STRUCTURES
struct Estimate
{
	ContigID cID;
	NumericID nID;
	int distance;
	int numPairs;
	
	friend std::ostream& operator<<(std::ostream& out, const Estimate& object)
	{
		out << object.cID << "," << object.distance << "," << object.numPairs;
		return out;
	} 
  
	friend std::istream& operator>> (std::istream& in, Estimate& object)
	{
		// Read 1 record from the stream
		std::string record;
		in >> record;
		
		// parse the record
		std::stringstream recss(record);
		std::stringstream convertor;
		std::string data;
	
		getline(recss, data, ',');
		convertor.str(data);
		convertor >> object.cID;
		
		object.nID = convertContigIDToNumericID(object.cID);
	
		getline(recss, data, ',');
		convertor.clear();
		convertor.str(data);	
		convertor >> object.distance;
		
		getline(recss, data, ',');
		convertor.clear();
		convertor.str(data);	
		convertor >> object.numPairs;
		return in;
	}  	
};

struct SimpleEdgeDesc
{
	ContigID contig;
	bool isRC;

	friend std::ostream& operator<<(std::ostream& out, const SimpleEdgeDesc& object)
	{
		out << object.contig << "," << object.isRC;
		return out;
	} 
  
	friend std::istream& operator>>(std::istream& in, SimpleEdgeDesc& object)
	{
		// Read 1 record from the stream
		std::string record;
		in >> record;
		
		// parse the record
		std::stringstream recss(record);
		std::stringstream convertor;
		std::string data;
	
		getline(recss, data, ',');
		convertor.str(data);
		convertor >> object.contig;
	
		getline(recss, data, ',');
		convertor.clear();
		convertor.str(data);	
		convertor >> object.isRC;
		return in;
	}    
};

typedef std::vector<Estimate> EstimateVector;

struct EstimateRecord
{
	NumericID refID;
	EstimateVector estimates[2];
};

// Estimate
void readEstimateRecord(std::ifstream& stream, EstimateRecord& er);

// Adjacency file
void loadGraphFromAdjFile(SimpleContigGraph* pGraph, ContigLengthMap& lengthMap, std::string file);
void parseAdjacencyLine(std::string& adjLine, ContigID currVert, SimpleContigGraph* pGraph);

// Length files
void loadContigLengths(std::string contigLenFile, ContigLengthMap& lengthMap);
int lookupLength(const ContigLengthMap& lengthMap, const ContigID& id);

// PDF loader
PDF loadPDF(std::string distCountFile, const int limit);

#endif
