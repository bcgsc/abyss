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

typedef std::vector<int> ContigLengthVec;
typedef DirectedGraph<SimpleContigData> SimpleContigGraph;

// STRUCTURES
struct Estimate
{
	LinearNumKey nID;
	int distance;
	int numPairs;
	int stdDev;
	
	friend std::ostream& operator<<(std::ostream& out, const Estimate& object)
	{
		out << object.nID << "," << object.distance << "," << object.numPairs << "," << object.stdDev;
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
		convertor >> object.nID;
	
		getline(recss, data, ',');
		convertor.clear();
		convertor.str(data);	
		convertor >> object.distance;
		
		getline(recss, data, ',');
		convertor.clear();
		convertor.str(data);	
		convertor >> object.numPairs;
		
		getline(recss, data, ',');
		convertor.clear();
		convertor.str(data);	
		convertor >> object.stdDev;		
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
	LinearNumKey refID;
	EstimateVector estimates[2];
};

// Estimate
void readEstimateRecord(std::ifstream& stream, EstimateRecord& er);

// Adjacency file
void loadGraphFromAdjFile(SimpleContigGraph* pGraph,  std::string& lengthFile, std::string adjFile);
void parseAdjacencyLine(std::string& adjLine, LinearNumKey currVert, SimpleContigGraph* pGraph);

// Length files
void loadContigLengths(std::string contigLenFile, ContigLengthVec& lengthMap);
int lookupLength(const ContigLengthVec& lengthMap, const LinearNumKey& id);

// PDF loader
PDF loadPDF(std::string distCountFile, const int limit);

// Convertor
LinearNumKey convertContigIDToLinearNumKey(const ContigID& id);



#endif
