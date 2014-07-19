#ifndef VECTOR_UTIL_H
#define VECTOR_UTIL_H_1

#include <vector>

template <typename T>
class make_vector {
private:
	std::vector<T> aVector;
public:
	typedef make_vector<T> mv;
	mv& operator<< (const T& val) {
		aVector.push_back(val);
		return *this;
	}
	operator std::vector<T>() const {
		return aVector;
	}
};

#endif
