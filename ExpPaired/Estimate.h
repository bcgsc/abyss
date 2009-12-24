#ifndef ESTIMATE_H
#define ESTIMATE_H 1

#include "PairUtils.h"
#include "StringUtil.h"
#include <cassert>
#include <cmath> // for ceilf
#include <iomanip>
#include <istream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <string>

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

typedef std::vector<Estimate> EstimateVector;

struct EstimateRecord
{
	LinearNumKey refID;
	EstimateVector estimates[2];

	/** Read the distance estimates for one contig. */
	friend std::istream& operator >>(std::istream& in,
			EstimateRecord& o)
	{
		o.estimates[SENSE].clear();
		o.estimates[ANTISENSE].clear();

		std::string id;
		in >> id;
		o.refID = convertContigIDToLinearNumKey(id);

		for (extDirection sense = SENSE;
				sense <= ANTISENSE; ++sense) {
			std::string s;
			std::getline(in, s, sense == SENSE ? ';' : '\n');
			std::istringstream ss(s);
			std::copy(std::istream_iterator<Estimate>(ss),
					std::istream_iterator<Estimate>(),
					std::back_inserter(o.estimates[sense]));
			assert(ss.eof());
		}

		return in;
	}

};

#endif
