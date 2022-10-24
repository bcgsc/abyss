/**
 * Functions for sequence manipulation.
 */
#ifndef BTLLIB_SEQ_HPP
#define BTLLIB_SEQ_HPP

#include "btllib/status.hpp"

#include <string>

namespace btllib {

extern const char COMPLEMENTS[256];
extern const char CAPITALS[256];

/**
 * Reverse complement a sequence in-place.
 *
 * @param seq Sequence to reverse complement.
 */
void
reverse_complement(std::string& seq);

/**
 * Obtain a reverse complement of the provided sequence. The argument sequence
 * is left untouched.
 *
 * @param seq Sequence to reverse complement.
 *
 * @return Reverse complemented sequence.
 */
std::string
get_reverse_complement(const std::string& seq);

} // namespace btllib

#endif