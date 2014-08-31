#ifndef VECTOR_UTIL_H
#define VECTOR_UTIL_H 1

#include <vector>

template <typename T>
class make_vector
{
private:

	std::vector<T> vec;

public:

	typedef make_vector<T> mv;
	typedef std::vector<T> sv;

	mv& operator<< (const T& val)
	{
		vec.push_back(val);
		return *this;
	}

	operator const sv&() const { return vec; }

	friend sv& operator+= (
			sv& vec1, const mv& vec2)
	{
		sv temp_vec = vec2;
		vec1.insert(vec1.end(), temp_vec.begin(), temp_vec.end());
		return vec1;
	}
};

#endif
