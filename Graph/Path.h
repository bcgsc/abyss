#ifndef PATH_H_
#define PATH_H_

#include <vector>
#include <string>
#include <sstream>
#include <climits>

enum PathSearchResult {
	FOUND_PATH = 0,
	TOO_MANY_PATHS,
	TOO_MANY_BRANCHES,
	EXCEEDED_MEM_LIMIT,
	NO_PATH
};

const char* PathSearchResultLabel[] = {
	"FOUND_PATH",
	"TOO_MANY_PATHS",
	"TOO_MANY_BRANCHES",
	"EXCEEDED_MEM_LIMIT",
	"NO_PATH"
};

enum Direction { FORWARD = 0, REVERSE };

const unsigned NO_LIMIT = UINT_MAX;

template <class Vertex> class Path : public std::vector<Vertex>
{
public:

	std::string str() {
		std::stringstream s;
		typename std::vector<Vertex>::iterator i = this->begin();
		for (; i != this->end(); i++) {
			if (i != this->begin())
				s << ",";
			s << *i;
		}
		return s.str();
	}
};

#endif
