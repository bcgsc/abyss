#ifndef ITERATOR_H
#define ITERATOR_H 1

#include <iterator>

/** A counting output iterator. */
class CountingOutputIterator
{
  public:
    explicit CountingOutputIterator(size_t& count) : m_count(count)
	{
		m_count = 0;
	}

	CountingOutputIterator& operator++()
	{
		++m_count;
		return *this;
	}

    CountingOutputIterator& operator*()
	{
		return *this;
	}

	template<typename T>
	CountingOutputIterator& operator=(const T&)
	{
		return *this;
	}

  private:
    size_t& m_count;
};

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
	typedef charT char_type;
	typedef traits traits_type;
	typedef std::basic_ostream<charT, traits> ostream_type;
	typedef void value_type;

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

/** Traits of an output iterator. */
template<class OutputIterator>
struct output_iterator_traits {
	typedef typename OutputIterator::value_type value_type;
};

/** Traits of a back_insert_iterator. */
template<class Container>
struct output_iterator_traits<
	std::back_insert_iterator<Container> >
{
	typedef typename Container::value_type value_type;
};

/** Traits of an ostream_iterator. */
template<class T, class charT, class traits>
struct output_iterator_traits<
	std::ostream_iterator<T, charT, traits> >
{
	typedef T value_type;
};

/** Traits of an affix_ostream_iterator. */
template<class T, class charT, class traits>
struct output_iterator_traits<
	affix_ostream_iterator<T, charT, traits> >
{
	typedef T value_type;
};

#endif
