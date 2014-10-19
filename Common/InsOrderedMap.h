/**
 * Partial implementation of a map class
 * read in multiple ways of ordering.
 *
 * index<0>: insertion order
 */

#ifndef INS_ORDERED_MAP_H
#define INS_ORDERED_MAP_H 1

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

using namespace boost::multi_index;

template <class F, class S>
class InsOrderedMap
{
private:

	struct aMap
	{
		F first;
		S second;
		aMap(const F first, const S second) : first(first), second(second) { }
	};

	typedef multi_index_container<
		aMap,
		indexed_by<
			random_access<>, // ra
			ordered_unique<member<aMap, F, &aMap::first> > // on
		>
	> aMap_cont;

	typedef typename aMap_cont::template nth_index<0>::type idx;
	typedef typename idx::iterator iit;

	aMap_cont ac;

	void insert_by_ra(const aMap_cont& other_ac)
	{
		const idx& i = other_ac.template get<0>();
		for (iit it = i.begin(), e=i.end(); it!=e; ++it)
			ac.push_back(aMap(it->first, it->second));
	}

public:

	InsOrderedMap() { }

	~InsOrderedMap() { }

	void push_back(const F& first, const S& second)
	{
		ac.push_back(aMap(first, second));
	}

	void insert(aMap_cont other_ac)
	{
		insert_by_ra(other_ac);
		other_ac.clear();
	}

	iit begin()
	{
		const idx& i = ac.template get<0>();
		return i.begin();
	}

	std::size_t size() { return ac.size(); }

	void erase(iit pair) { ac.erase(pair); }

	void clear() { ac.clear(); }

	bool empty() { return ac.empty(); }

	const F& getFirst(const iit& pair) { return pair->first; }

	const S& getSecond(const iit& pair) { return pair->second; }

	const aMap_cont& getAC() { return ac; }
};

#endif
