#ifndef AFFIXITERATOR_H
#define AFFIXITERATOR_H 1

#include <iterator>

/** An output stream iterator, like ostream_iterator, that outputs the
 * a delimiter before and after the item.
 */
template <class T,
		class charT = char,
		class traits = std::char_traits<charT> >

class affix_ostream_iterator :
    public std::iterator<std::output_iterator_tag,
		void, void, void, void>
{
public:
	/** The type of element inserted into the stream. Note that the
	 * standard ostream_iterator does not have this typedef. */
	typedef T value_type;

	typedef charT char_type;
	typedef traits traits_type;
	typedef std::basic_ostream<charT, traits> ostream_type;

	affix_ostream_iterator(ostream_type& s,
			charT const *prefix, charT const *suffix = NULL)
		: os(&s), prefix(prefix), suffix(suffix)
	{}

	affix_ostream_iterator<T, charT,
		traits>& operator =(T const &item)
	{
		if (prefix != NULL)
			*os << prefix;
		*os << item;
		if (suffix != NULL)
			*os << suffix;
		return *this;
	}

	affix_ostream_iterator<T, charT, traits> &operator*() {
		return *this;
	}
	affix_ostream_iterator<T, charT, traits> &operator++() {
		return *this;
	}
	affix_ostream_iterator<T, charT, traits> &operator++(int) {
		return *this;
	}

private:
	std::basic_ostream<charT, traits> *os;
	charT const* prefix;
	charT const* suffix;
};

#endif
