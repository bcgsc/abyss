#ifndef PREFIXITERATOR_H
#define PREFIXITERATOR_H 1

#include <iterator>

/** An output stream iterator, like ostream_iterator, that outputs the
 * delimiter before the item rather than after.
 */
template <class T,
		class charT = char,
		class traits = std::char_traits<charT> >

class prefix_ostream_iterator :
    public std::iterator<std::output_iterator_tag,
		void, void, void, void>
{
public:
	typedef charT char_type;
	typedef traits traits_type;
	typedef std::basic_ostream<charT, traits> ostream_type;

	prefix_ostream_iterator(ostream_type& s)
		: os(&s), delimiter(0)
	{}

	prefix_ostream_iterator(ostream_type& s, charT const *d)
		: os(&s), delimiter(d)
	{}

	prefix_ostream_iterator<T, charT,
		traits>& operator =(T const &item)
	{
		// Here's the only real change from ostream_iterator:
		// Normally, the '*os << item;' would come before the 'if'.
		if (delimiter != 0)
			*os << delimiter;
		*os << item;
		return *this;
	}

	prefix_ostream_iterator<T, charT, traits> &operator*() {
		return *this;
	}
	prefix_ostream_iterator<T, charT, traits> &operator++() {
		return *this;
	}
	prefix_ostream_iterator<T, charT, traits> &operator++(int) {
		return *this;
	}

private:
	std::basic_ostream<charT, traits> *os;
	charT const* delimiter;
};

#endif
