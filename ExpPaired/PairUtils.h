#ifndef PAIRUTILS_H
#define PAIRUTILS_H 1

#include "Dictionary.h"
#include "StringUtil.h"
#include <cassert>
#include <cmath> // for ceilf
#include <iomanip>
#include <iostream>
#include <stdint.h>
#include <vector>

typedef std::vector<int> ContigLengthVec;

typedef uint32_t LinearNumKey;

extern Dictionary g_contigIDs;

/** Distance estimate between two contigs. */
struct Estimate
{
	LinearNumKey nID;
	int distance;
	unsigned numPairs;
	float stdDev;
	bool isRC;
	
	friend std::ostream& operator<<(std::ostream& out,
			const Estimate& o)
	{
		return out << g_contigIDs.key(o.nID)
			<< (o.isRC ? '-' : '+') << ","
			<< o.distance << ","
			<< o.numPairs << ","
			<< std::fixed << std::setprecision(1) << o.stdDev;
	}

	friend std::istream& operator>>(std::istream& in, Estimate& o)
	{
		in >> std::ws;
		std::string sID;
		char comma0, comma1;
		getline(in, sID, ',')
			>> o.distance >> comma0
			>> o.numPairs >> comma1
			>> o.stdDev;
		if (in) {
			assert(comma0 == ',' && comma1 == ',');
			char c = chop(sID);
			assert(c == '+' || c == '-');
			o.isRC = c == '-';
			o.nID = g_contigIDs.serial(sID);
		}
		return in;
	}
};

/** Return the allowed error for the given estimate. */
static inline unsigned allowedError(float stddev)
{
	/** The number of standard deviations. */
	const int NUM_SIGMA = 3;

	/**
	 * Additional constant error. The error expected that does not
	 * vary with the number of samples.
	 */
	const unsigned CONSTANT_ERROR = 6;

	return (unsigned)ceilf(NUM_SIGMA * stddev + CONSTANT_ERROR);
}


typedef std::string ContigID;

struct SimpleEdgeDesc
{
	ContigID contig;
	bool isRC;

	SimpleEdgeDesc() { }
	SimpleEdgeDesc(ContigID contig, bool isRC)
		: contig(contig), isRC(isRC) { }

	friend std::ostream& operator<<(std::ostream& out,
			const SimpleEdgeDesc& o)
	{
		return out << o.contig << (o.isRC ? '-' : '+');
	}

	friend std::istream& operator>>(std::istream& in,
			SimpleEdgeDesc& o)
	{
		if (in >> o.contig) {
			char c = chop(o.contig);
			assert(c == '+' || c == '-');
			o.isRC = c == '-';
		}
		return in;
	}
};

typedef std::vector<Estimate> EstimateVector;

struct EstimateRecord;

struct EstimateRecord
{
	LinearNumKey refID;
	EstimateVector estimates[2];
	friend std::istream& operator >>(std::istream& in,
			EstimateRecord& er);
};

void loadContigLengths(const std::string& path,
		ContigLengthVec& lengths);

LinearNumKey convertContigIDToLinearNumKey(const ContigID& id);

#endif
