#ifndef PATH_H_
#define PATH_H_

#include <vector>
#include <string>
#include <sstream>

enum PathSearchResult {
	FOUND_PATH,
	TOO_MANY_PATHS,
	TOO_MANY_BRANCHES,
	NO_PATH
};

enum Direction { FORWARD, REVERSE };

const unsigned NO_LIMIT = 0;

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

/*
template <class Graph> class PathList :
	public std::vector<std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> > {};
*/

#endif
