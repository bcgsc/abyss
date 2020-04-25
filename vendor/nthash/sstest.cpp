#include <iostream>
#include <string>
#include <vector>
#include "ssHashIterator.hpp"

using namespace std;

int main()
{
    /* test sequence */
    std::string seq = "GAGTGTCAAACATTCAGACAACAGCAGGGGTGCTCTGGAATCCTATGTGAGGAACAAACATTCAGGCCACAGTAG";

    /* seed is the spaced seed */
    std::string seedStr("110001100111000001110100110110100010011101");

    std::vector<bool> seed;
    for(auto p : seedStr)
        seed.push_back(p=='1');

    ssHashIterator ssitr(seq, seed, seed.size());

    while (ssitr != ssitr.end()) {
        std::cout << *ssitr << std::endl;
        ++ssitr;
    }

    return 0;
}