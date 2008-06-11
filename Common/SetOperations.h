#ifndef SETOPERATIONS_H
#define SETOPERATIONS_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <iterator>

#include <set>

namespace SetOperations
{

template<typename K>
std::ostream& printSet(const std::set<K>& s)
{
	std::cout << "{ ";
	std::copy(s.begin(), s.end(),
			std::ostream_iterator<K>(std::cout, " "));
	std::cout << "}";
	return std::cout;
}

template<typename K>
std::ostream& printPath(const std::vector<K>& s)
{
	std::cout << "[ ";
	std::copy(s.begin(), s.end(),
			std::ostream_iterator<K>(std::cout, " "));
	std::cout << "]";
	return std::cout;
}


template<typename K>
void difference(const std::set<K>& s1, const std::set<K>& s2, std::set<K>& result)
{ 
	std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(result, result.begin()));
}

template<typename K>
void removeElements(std::set<K>& s1, const std::set<K>& s2)
{
	std::set<K> result;
	difference(s1, s2, result);
	s1 = result;
}

};

#endif
