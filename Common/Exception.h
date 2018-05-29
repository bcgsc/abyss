#ifndef _EXCEPTION_H_
#define _EXCEPTION_H_ 1

/* `noexcept` is only recognized in C++11 or later. See https://stackoverflow.com/questions/24567173/backwards-compatible-noexceptfalse-for-destructors */
#if __cplusplus >= 201103L
#define NOEXCEPT noexcept
#else
#define NOEXCEPT
#endif

#endif
