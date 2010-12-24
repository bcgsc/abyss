#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H 1

#include <functional>

/** A functor that always returns true. */
template <typename Arg>
struct True : std::unary_function<Arg, bool> {
	bool operator()(const Arg&) const { return true; }
};

/** A functor to adapt a pointer to a member variable. */
template <typename Arg, typename Result>
struct MemVar : std::unary_function<Arg, Result>
{
	MemVar(Result Arg::* p) : m_p(p) { }
	Result operator()(const Arg& arg) const { return arg.*m_p; }
  private:
	Result Arg::* m_p;
};

/** Return a functor to adapt a pointer to a member variable. */
template <typename Arg, typename Result>
MemVar<Arg, Result> mem_var(Result Arg::* p)
{
	return MemVar<Arg, Result>(p);
}

/** A functor, f1(f2(x)). */
template <typename F1, typename F2>
struct unary_compose : std::unary_function <
		typename F2::argument_type, typename F1::result_type>
{
	unary_compose(const F1& f1, const F2& f2) : f1(f1), f2(f2) { }
	typename F1::result_type operator()(
			const typename F2::argument_type& x) const
	{
		return f1(f2(x));
	}
  private:
	F1 f1;
	F2 f2;
};

/** Return a functor, f1(f2(x)). */
template <typename F1, typename F2>
unary_compose<F1, F2> compose1(const F1& f1, const F2& f2)
{
	return unary_compose<F1, F2>(f1, f2);
}

/** A functor, f(g1(x), g2(x)). */
template <typename F, typename G1, typename G2>
struct binary_compose : std::unary_function <
		typename G1::argument_type, typename F::result_type>
{
	binary_compose(const F& f, const G1& g1, const G2& g2)
		: f(f), g1(g1), g2(g2) { }
	typename G1::result_type operator()(
			const typename G1::argument_type& x) const
	{
		return f(g1(x), g2(x));
	}
  private:
	F f;
	G1 g1;
	G2 g2;
};

/** Return a functor, f(g1(x), g2(x)). */
template <typename F, typename G1, typename G2>
binary_compose<F, G1, G2> compose2(
		const F& f, const G1& g1, const G2& g2)
{
	return binary_compose<F, G1, G2>(f, g1, g2);
}

#endif
