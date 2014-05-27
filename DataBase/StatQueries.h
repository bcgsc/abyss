/**
 * Create tables and insert values into database repository.
 */

#ifndef STATQUERIES_H
#define STATQUERIES_H 1

#include "DB.h"

template <typename D>
static inline void add2db (D& db, const std::string& key, const float& value)
{
	db.statMap[key] = value;
}

template <typename D, typename Map>
static inline void add2db (D& db, const Map& m)
{
	db.statMap.insert (m.begin(), m.end());
}

#endif
