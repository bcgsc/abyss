#ifndef PATH_H_
#define PATH_H_

#include <string>
#include <sstream>
#include <climits>
#include <deque>
#include <cassert>

enum PathSearchResult {
	FOUND_PATH = 0,
	TOO_MANY_PATHS,
	TOO_MANY_BRANCHES,
	PATH_CONTAINS_CYCLE,
	MAX_COST_EXCEEDED,
	EXCEEDED_MEM_LIMIT,
	NO_PATH
};

const char* PathSearchResultLabel[] = {
	"FOUND_PATH",
	"TOO_MANY_PATHS",
	"TOO_MANY_BRANCHES",
	"PATH_CONTAINS_CYCLE",
	"MAX_COST_EXCEEDED",
	"EXCEEDED_MEM_LIMIT",
	"NO_PATH"
};

enum Direction { FORWARD = 0, REVERSE };

inline static const char* directionStr(Direction dir)
{
	switch(dir) {
	case FORWARD:
		return "FORWARD";
	case REVERSE:
		return "REVERSE";
	default:
		assert(false);
	}
}

const unsigned NO_LIMIT = UINT_MAX;

template <class Vertex> class Path : public std::deque<Vertex>
{
public:

	std::string str() {
		std::stringstream s;
		typename std::deque<Vertex>::iterator i = this->begin();
		for (; i != this->end(); i++) {
			if (i != this->begin())
				s << ",";
			s << *i;
		}
		return s.str();
	}
};

#endif
