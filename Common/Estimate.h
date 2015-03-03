#ifndef ESTIMATE_H
#define ESTIMATE_H 1

#include "Common/ContigProperties.h" // for Distance
#include "ContigID.h"
#include "ContigNode.h"
#include "Graph/Options.h" // for opt::k
#include "IOUtil.h"
#include <cassert>
#include <cmath> // for ceilf
#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace opt {
	/** The acceptable error of a distance estimate. */
	extern unsigned distanceError;
}

/** An estimate of the distance between two contigs. */
struct DistanceEst
{
	int distance;
	unsigned numPairs;
	float stdDev;

	DistanceEst() : distance(-opt::k + 1), numPairs(0), stdDev(0) { }

	DistanceEst(int distance)
		: distance(distance), numPairs(0), stdDev(distance) { }

	DistanceEst(int distance, unsigned numPairs, float stdDev)
		: distance(distance), numPairs(numPairs), stdDev(stdDev) { }

	bool operator==(const DistanceEst& o) const
	{
		return distance == o.distance
			&& numPairs == o.numPairs
			&& stdDev == o.stdDev;
	}

	friend std::ostream& operator<<(std::ostream& out,
			const DistanceEst& o)
	{
		if (opt::format != DIST) {
			out << "d=" << o.distance;
			if (o.stdDev > 0 || o.numPairs > 0)
				out << " e=" << std::fixed << std::setprecision(1)
					<< o.stdDev
					<< " n=" << o.numPairs;
			return out;
		} else
			return out << o.distance << ',' << o.numPairs << ','
				<< std::fixed << std::setprecision(1) << o.stdDev;
	}

	friend std::istream& operator>>(std::istream& in, DistanceEst& o)
	{
		if (in >> std::ws && in.peek() == 'd') {
			if (!(in >> expect("d =") >> o.distance >> std::ws))
				return in;
			if (in.peek() == ']') {
				o.stdDev = o.numPairs = 0;
				return in;
			} else if (in.peek() == ',')
				return in >> expect(", e =") >> o.stdDev
					>> expect(", n =") >> o.numPairs;
			else
				return in >> expect(" e =") >> o.stdDev
					>> expect(" n =") >> o.numPairs;
		} else
			return in >> o.distance >> expect(",")
				>> o.numPairs >> expect(",")
				>> o.stdDev;
	}
};

/** Return the better of two distance estimates.
 * Return the estimate whose error is least, or when the errors are
 * equal, return the larger distance estimate.
 * Add the number of pairs.
 */
struct BetterDistanceEst
{
	NoProperty operator()(
			const NoProperty& a, const NoProperty&) const
	{
		return a;
	}

	Distance operator()(
			const Distance& a, const Distance& b) const
	{
		return a.distance > b.distance ? a : b;
	}

	DistanceEst operator()(
			const DistanceEst& a, const DistanceEst& b) const
	{
		bool which = a.stdDev != b.stdDev ? a.stdDev < b.stdDev
			: a.distance > b.distance;
		DistanceEst x = which ? a : b;
		x.numPairs = a.numPairs + b.numPairs;
		return x;
	}
};

/** Merge two distance estimates. */
struct MergeDistanceEst
{
	DistanceEst operator()(
			const DistanceEst& a, const DistanceEst& b) const
	{
		int x1 = a.distance, x2 = b.distance;
		double v1 = a.stdDev * a.stdDev, v2 = b.stdDev * b.stdDev;
		DistanceEst x;
		x.distance = (int)round((x1 * v2 + x2 * v1) / (v1 + v2));
		x.stdDev = sqrt(v1 * v2 / (v1 + v2));
		x.numPairs = a.numPairs + b.numPairs;
		return x;
	}
};

/** Return the allowed error for the given estimate. */
static inline unsigned allowedError(float stddev)
{
	/** The number of standard deviations. */
	const int NUM_SIGMA = 3;
	return (unsigned)ceilf(NUM_SIGMA * stddev + opt::distanceError);
}

/** Distance estimates to and from a particular contig. */
struct EstimateRecord
{
	typedef std::pair<ContigNode, DistanceEst> Estimate;
	typedef std::vector<Estimate> Estimates;

	ContigID refID;
	Estimates estimates[2];

	/** Read the distance estimates for one contig. */
	friend std::istream& operator >>(std::istream& in,
			EstimateRecord& o)
	{
		o.estimates[false].clear();
		o.estimates[true].clear();

		std::string name;
		in >> name;
		if (!in)
			return in;
		o.refID = ContigID(get(g_contigNames, name));

		for (int rc = false; rc <= true; ++rc) {
			std::string s;
			std::getline(in, s, !rc ? ';' : '\n');
			std::istringstream ss(s);
			for (Estimate ep; getline(ss >> std::ws, s, ',');) {
				ep.first = find_vertex(s, g_contigNames);
				if (ss >> ep.second)
					o.estimates[rc].push_back(ep);
			}
			assert(ss.eof());
		}

		return in;
	}

};

namespace std {
	template<>
	inline void swap(EstimateRecord&, EstimateRecord&)
	{
		assert(false);
	}
}

#endif
