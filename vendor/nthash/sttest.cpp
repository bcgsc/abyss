#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "stHashIterator.hpp"
using namespace std;

int main(int argc, const char *argv[])
{
    size_t countOnes=0;
    
    std::vector<std::string> seedString;
    seedString.push_back("110001100111000001110100110110100010011101");
    seedString.push_back("000001111100101110111000001011010100011110");
    seedString.push_back("011110001010110100000111011101001111100000");
    seedString.push_back("101110010001011011001011100000111001100011");
    
    
    /* test sequence */
    //std::string seq = "GAGTGTCAAACATTCAGACAACAGCAGGGGTGCTCTGGAATCCTATGTGAGGAACAAACATTCAGGCCACAGTAG";
    
    /* h is the number of spaced seeds and k is the length of spaced seeds */
    unsigned h, k;
    h = seedString.size();
    k = seedString[0].size();
    
    std::vector<std::vector<unsigned> > seedSet = parseSeed(seedString);
    
    clock_t sTime = clock();
    std::ifstream in(argv[argc-1]);
    bool good = true;
    for(string seq, hseq; good;) {
        good = static_cast<bool>(getline(in, hseq));
        good = static_cast<bool>(getline(in, seq));
        good = static_cast<bool>(getline(in, hseq));
        good = static_cast<bool>(getline(in, hseq));
        if(good){
            stHashIterator ssitr(seq, seedSet, h, k);
            while (ssitr != ssitr.end()) {
                if((ssitr.strandArray())[0])
                    ++countOnes;
                ++ssitr;
            }
        }
        
    }

    std::cerr<< countOnes << std::endl;
    std::cerr << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
    
    
    return 0;
}
