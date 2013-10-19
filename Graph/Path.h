#ifndef PATH_H_
#define PATH_H_

enum PathSearchResult {
	FOUND_PATH,
	TOO_MANY_PATHS,
	TOO_MANY_BRANCHES,
	NO_PATH
};

enum Direction { FORWARD, REVERSE };

const unsigned NO_LIMIT = 0;

#endif
