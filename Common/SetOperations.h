#ifndef SETOPERATIONS_H
#define SETOPERATIONS_H

#include <ostream>
#include <iterator>
#include <map>
#include <vector>

namespace SetOperations
{

template<typename K, typename D>
std::ostream& printMap(const std::map<K,D>& s)
{
	std::cout << "{ ";
	for(typename std::map<K,D>::const_iterator iter = s.begin(); iter != s.end(); ++iter)
	{
		std::cout << iter->first << "," << iter->second << " ";
	}

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

}

#endif
