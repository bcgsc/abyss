#ifndef ESTIMATE_H
#define ESTIMATE_H 1

#include "ContigID.h"
#include "ContigNode.h"
#include <cassert>
#include <cmath> // for ceilf
#include <iomanip>
#include <istream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <string>

namespace opt {
	extern int dot;
}

/** Distance estimate between two contigs. */
struct Estimate
{
	ContigNode contig;
	int distance;
	unsigned numPairs;
	float stdDev;

	friend std::ostream& operator<<(std::ostream& out,
			const Estimate& o)
	{
		if (opt::dot)
			return out << '"' << o.contig << "\" ["
				"d=" << o.distance << " "
				"e=" << std::fixed << std::setprecision(1)
					<< o.stdDev << " "
				"n=" << o.numPairs << ']';
		else
			return out << o.contig << ','
				<< o.distance << ','
				<< o.numPairs << ','
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
		if (in)
			o.contig = ContigNode(sID);
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
	ContigID refID;
	EstimateVector estimates[2];

	/** Read the distance estimates for one contig. */
	friend std::istream& operator >>(std::istream& in,
			EstimateRecord& o)
	{
		o.estimates[false].clear();
		o.estimates[true].clear();

		ContigID id;
		in >> id;
		if (!in)
			return in;
		o.refID = id;

		for (int rc = false; rc <= true; ++rc) {
			std::string s;
			std::getline(in, s, !rc ? ';' : '\n');
			std::istringstream ss(s);
			std::copy(std::istream_iterator<Estimate>(ss),
					std::istream_iterator<Estimate>(),
					std::back_inserter(o.estimates[rc]));
			assert(ss.eof());
		}

		return in;
	}

};

#endif
