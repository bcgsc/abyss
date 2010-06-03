#ifndef CONTIGLENGTH_H
#define CONTIGLENGTH_H 1

#include <istream>
#include <string>
#include <vector>

std::vector<unsigned> readContigLengths(std::istream& in);
std::vector<unsigned> readContigLengths(const std::string& path);

#endif
