#ifndef SENSE_H
#define SENSE_H 1

#include <cassert>

enum extDirection
{
	SENSE = 0,
	ANTISENSE = 1,
	NUM_DIRECTIONS
};

static inline extDirection operator !(extDirection dir)
{
	return dir == SENSE ? ANTISENSE : SENSE;
}

static inline extDirection& operator ++(extDirection& dir)
{
	assert(dir == SENSE || dir == ANTISENSE);
	return dir = extDirection(dir + 1);
}

#endif
