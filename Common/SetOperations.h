#ifndef SETOPERATIONS_H
#define SETOPERATIONS_H

#include <set>
#include <iostream>

namespace SetOperations
{

template<typename K>
void printSet(const std::set<K>& s)
{
	std::cout << "[ ";
	for(typename std::set<K>::const_iterator iter = s.begin(); iter != s.end(); ++iter)
	{
		std::cout << *iter << " ";
	}
	std::cout << "]";
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
