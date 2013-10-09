#ifndef PATH_H_
#define PATH_H_

enum PathSearchResult { 
	FOUND_PATH, 
	TOO_MANY_PATHS, 
	TOO_MANY_BRANCHES, 
	NO_PATH 
};

extern const int NO_LIMIT = -1;

#endif
