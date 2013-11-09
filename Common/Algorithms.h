#ifndef ALGORITHMS_H
#define ALGORITHMS_H 1

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

/** Apply function f to each element in the range [first, last) for
 * which the predicate p is true.
 */
template<class It, class Function, class Pred>
Function for_each_if(It first, It last, Function f, Pred p)
{
	for (; first != last; ++first)
		if (p(*first))
			f(*first);
	return f;
}

/** Copies each element in the range [first, last) into the range
 * starting at result for which the predicate p is true. */
template<class InputIt, class OutputIt, class Pred>
OutputIt copy_if(InputIt first, InputIt last, OutputIt result,
		Pred pred)
{
	for (; first != last; ++first) {
		if (pred(*first)) {
			*result = *first;
			++result;
		}
	}
	return result;
}

/** Sorts the elements in the range [first,last) ordered by the value
 * returned by the unary function op, which is called once for each
 * element in the range [first,last). The copy constructor of the
 * value_type of It is not used.
 */
template <class It, class Op>
void sort_by_transform(It first, It last, Op op)
{
	typedef typename std::iterator_traits<It>::difference_type
		size_type;
	typedef typename Op::result_type key_type;

	size_type n = last - first;
	std::vector< std::pair<key_type, size_type> > keys;
	keys.reserve(n);
	for (It it = first; it != last; ++it)
		keys.push_back(std::make_pair(op(*it), it - first));
	sort(keys.begin(), keys.end());

	for (size_type i = 0; i < n; i++) {
		size_type j = keys[i].second;
		while (j < i)
			j = keys[j].second;
		if (i != j)
			std::swap(first[i], first[j]);
		keys[i].second = j;
	}
}

#endif
