#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H 1

#include <functional>

/** Always return true. */
template <typename Arg>
struct True : std::unary_function<Arg, bool> {
	bool operator()(const Arg&) const { return true; }
};

/** Compose two unary functions. */
template <typename F1, typename F2>
struct unary_compose : std::unary_function <
        typename F2::argument_type, typename F1::result_type>
{
    unary_compose(F1 f1, F2 f2) : f1(f1), f2(f2) {}

    typename F1::result_type operator()(typename F2::argument_type x)
    {
		return f1(f2(x));
	}

private:
    F1 f1;
    F2 f2;
};

/** Compose two unary functions. */
template <typename F1, typename F2>
unary_compose<F1, F2> compose1(F1 f1, F2 f2)
{
	return unary_compose<F1, F2>(f1, f2);
}

/** Return the logical and of two unary predicates. */
template <typename F1, typename F2>
struct PredicateAnd
	: std::unary_function<typename F1::argument_type, bool>
{
	PredicateAnd(const F1& f1, const F2& f2) : f1(f1), f2(f2) { }
	bool operator()(const typename F1::argument_type& x) const
	{
		return f1(x) && f2(x);
	}
  private:
	F1 f1;
	F2 f2;
};

#endif
