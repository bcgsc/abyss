#ifndef CONTIGPROPERTIES_H
#define CONTIGPROPERTIES_H 1

#include <istream>
#include <ostream>

/** Contig properties. */
struct ContigProperties {
	unsigned length;
	unsigned coverage;

	ContigProperties() { }
	ContigProperties(unsigned length, unsigned coverage)
		: length(length), coverage(coverage) { }

	friend std::ostream& operator <<(std::ostream& out,
			const ContigProperties& o)
	{
		return out << ' ' << o.length << ' ' << o.coverage;
	}

	friend std::istream& operator >>(std::istream& in,
			ContigProperties& o)
	{
		return in >> o.length >> o.coverage;
	}
};

#endif
