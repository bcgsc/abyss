/*
 *
 * nttest.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */


#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <unistd.h>
#include <getopt.h>
#include "seqgen.hpp"
#include "BloomFilter.hpp"
#include "nthash.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "nttest"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamid Mohamadi.\n"
    "Copyright 2016 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... QUERY\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";

namespace opt {
unsigned threads=1;
unsigned kmerLen=50;
unsigned ibits = 8;
unsigned nhash=1;
unsigned nz;
size_t nquery;
size_t squery;
size_t ngene;
size_t sgene;
unsigned method;
unsigned maxitr;
bool fastq = false;
bool inpFlag = false;
bool uniformity = false;
}

using namespace std;

static const char shortopts[] = "k:b:h:j:q:l:t:g:m:a:iu";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 'j' },
    { "kmer",	required_argument, NULL, 'k' },
    { "qnum",	required_argument, NULL, 'q' },
    { "qlen",	required_argument, NULL, 'l' },
    { "bit",	required_argument, NULL, 'b' },
    { "hash",	required_argument, NULL, 'h' },
    { "tnum",	required_argument, NULL, 't' },
    { "tlen",	required_argument, NULL, 'g' },
    { "max",	required_argument, NULL, 'm' },
    { "alg",	required_argument, NULL, 'a' },
    { "input",	no_argument, NULL, 'i' },
    { "uniformity",	no_argument, NULL, 'u' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

static const string itm[]= {"ntbase","nthash","city"};

void getFtype(const char *fName) {
    std::ifstream in(fName);
    std::string hLine;
    bool good=static_cast<bool>(getline(in,hLine));
    in.close();
    if(!good) {
        std::cerr<<"Error in reading file: "<<fName<<"\n";
        exit(EXIT_FAILURE);
    }
    if(hLine[0]=='>')
        opt::fastq=false;
    else if (hLine[0]=='@')
        opt::fastq=true;
    else {
        std::cerr<<"Error in file format: "<<fName<<"\n";
        exit(EXIT_FAILURE);
    }
}

bool getSeq(std::ifstream &uFile, std::string &line) {
    bool good=false;
    std::string hline;
    line.clear();
    if(opt::fastq) {
        good=static_cast<bool>(getline(uFile, hline));
        good=static_cast<bool>(getline(uFile, line));
        good=static_cast<bool>(getline(uFile, hline));
        good=static_cast<bool>(getline(uFile, hline));
    }
    else {
        do {
            good=static_cast<bool>(getline(uFile, hline));
            if(hline[0]=='>'&&!line.empty()) break;// !line.empty() for the first rec
            if(hline[0]!='>')line+=hline;
        } while(good);
        if(!good&&!line.empty())
            good=true;
    }
    return good;
}

static const unsigned char b2r[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //0
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //1
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //2
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //3
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //4   'A' 'C' 'G'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //5   'T'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //6   'a' 'c' 'g'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //7   't'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //8
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //9
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //10
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //11
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //12
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //13
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //14
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //15
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

void getCanon(std::string &bMer) {
    int p=0, hLen=(opt::kmerLen-1)/2;
    while (bMer[p] == b2r[(unsigned char)bMer[opt::kmerLen-1-p]]) {
        ++p;
        if(p>=hLen) break;
    }
    if (bMer[p] > b2r[(unsigned char)bMer[opt::kmerLen-1-p]]) {
        for (int lIndex = p, rIndex = opt::kmerLen-1-p; lIndex<=rIndex; ++lIndex,--rIndex) {
            char tmp = b2r[(unsigned char)bMer[rIndex]];
            bMer[rIndex] = b2r[(unsigned char)bMer[lIndex]];
            bMer[lIndex] = tmp;
        }
    }
}

void loadSeq(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        myFilter.insert(seq.c_str()+i);
    }
}

void loadSeqr(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    uint64_t fhVal, rhVal;
    myFilter.insert(seq.c_str(), fhVal, rhVal);
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        myFilter.insert(fhVal, rhVal, seq[i-1], seq[i+opt::kmerLen-1]);
    }
}

void loadSeqm(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        myFilter.insertMur(kmer.c_str());
    }
}

void loadSeqc(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        myFilter.insertCit(kmer.c_str());
    }
}

void loadSeqx(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        myFilter.insertXxh(kmer.c_str());
    }
}

void loadBf(BloomFilter &myFilter, const char* faqFile) {
    getFtype(faqFile);
    ifstream uFile(faqFile);
    bool good = true;
    
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    
    for(string line; good;) {

        #ifdef _OPENMP
        #pragma omp critical(uFile)
        #endif

        good = getSeq(uFile, line);
        if(good) {
            if(itm[opt::method]=="city")
                loadSeqc(myFilter, line);
            else if(itm[opt::method]=="murmur")
                loadSeqm(myFilter, line);
            else if(itm[opt::method]=="xxhash")
                loadSeqx(myFilter, line);
            else if(itm[opt::method]=="ntbase")
                loadSeq(myFilter, line);
            else if(itm[opt::method]=="nthash")
                loadSeqr(myFilter, line);
        }
    }
    uFile.close();
}

void querySeqr(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
    uint64_t fhVal, rhVal;
    if(myFilter.contains(seq.c_str(), fhVal, rhVal)) {
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        ++fHit;
    }
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        if(myFilter.contains(fhVal, rhVal, seq[i-1], seq[i+opt::kmerLen-1])) {
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            ++fHit;
        }
    }
}

void querySeq(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        if(myFilter.contains(seq.c_str()+i)) {
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            ++fHit;
        }
    }
}

void querySeqm(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        if(myFilter.containsMur(kmer.c_str())) {
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            ++fHit;
        }
    }
}

void querySeqc(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        if(myFilter.containsCit(kmer.c_str())) {
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            ++fHit;
        }
    }
}

void querySeqx(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        if(myFilter.containsXxh(kmer.c_str())) {
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            ++fHit;
        }
    }
}

void queryBf(BloomFilter &myFilter, const char* faqFile) {
    getFtype(faqFile);
    ifstream uFile(faqFile);
    size_t fHit=0,totKmer=0;
    bool good = true;
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for(string line; good;) {
        #ifdef _OPENMP
        #pragma omp critical(uFile)
        #endif
        good = getSeq(uFile, line);
        if(good) {
            if(itm[opt::method]=="city")
                querySeqc(myFilter, line, fHit);
            else if(itm[opt::method]=="murmur")
                querySeqm(myFilter, line, fHit);
            else if(itm[opt::method]=="xxhash")
                querySeqx(myFilter, line, fHit);
            else if(itm[opt::method]=="ntbase")
                querySeq(myFilter, line, fHit);
            else if(itm[opt::method]=="nthash")
                querySeqr(myFilter, line, fHit);
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            totKmer+=opt::squery-opt::kmerLen+1;
        }
    }
    uFile.close();
    cerr << "tkmer=" << totKmer << " ";
    cerr << "fhits=" << fHit << " %" << setprecision(4) << fixed << (double)fHit/(double)totKmer << " ";
}

void hashSeqb(const string & seq) {
    for (size_t i = 0; i < seq.length() - opt::kmerLen + 1; i++) {
        if(NTC64(seq.c_str()+i, opt::kmerLen)) opt::nz++;
    }
}

void hashSeqr(const string & seq) {
    uint64_t fhVal,rhVal,hVal;
    hVal = NTC64(seq.c_str(), opt::kmerLen, fhVal, rhVal);
    if(hVal)opt::nz++;
    for (size_t i = 1; i < seq.length() - opt::kmerLen + 1; i++) {
        hVal = NTC64(seq[i-1], seq[i-1+opt::kmerLen], opt::kmerLen, fhVal, rhVal);
        if(hVal)opt::nz++;
    }
}

void hashSeqx(const string & seq) {
    for (size_t i = 0; i < seq.length() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        if(XXH64(kmer.c_str(), opt::kmerLen, 0))++opt::nz;
    }
}

void hashSeqc(const string & seq) {
    for (size_t i = 0; i < seq.length() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        if (CityHash64(kmer.c_str(), opt::kmerLen))++opt::nz;
    }
}

void hashSeqm(const string & seq) {
    for (size_t i = 0; i < seq.length() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        if(MurmurHash64A(kmer.c_str(), opt::kmerLen, 0))++opt::nz;
    }
}

void hashSeqbM(const string & seq) {
    uint64_t hVal[opt::nhash];
    for (size_t i = 0; i < seq.length() - opt::kmerLen + 1; i++) {
        NTMC64(seq.c_str()+i, opt::kmerLen, opt::nhash, hVal);
        for(unsigned i=0; i<opt::nhash; i++)
            if(hVal[i])opt::nz++;
    }
}

void hashSeqrM(const string & seq) {
    uint64_t hVal[opt::nhash], fhVal, rhVal;
    NTMC64(seq.c_str(), opt::kmerLen, opt::nhash, fhVal, rhVal, hVal);
    for(unsigned h=0; h<opt::nhash; h++) if(hVal[h])opt::nz++;
    for (size_t i = 1; i < seq.length() - opt::kmerLen + 1; i++) {
        NTMC64(seq[i-1], seq[i-1+opt::kmerLen], opt::kmerLen, opt::nhash, fhVal, rhVal, hVal);
        for(unsigned h=0; h<opt::nhash; h++) if(hVal[h])opt::nz++;
    }
}

void hashSeqxM(const string & seq) {
    for (size_t i = 0; i < seq.length() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        for(unsigned h=0; h<opt::nhash; h++)
            if(XXH64(kmer.c_str(), opt::kmerLen, h))++opt::nz;
    }
}

void hashSeqcM(const string & seq) {
    for (size_t i = 0; i < seq.length() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        for(unsigned h=0; h<opt::nhash; h++)
            if (CityHash64WithSeed(kmer.c_str(), opt::kmerLen, h))++opt::nz;
    }
}

void hashSeqmM(const string & seq) {
    for (size_t i = 0; i < seq.length() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        for(unsigned h=0; h<opt::nhash; h++)
            if(MurmurHash64A(kmer.c_str(), opt::kmerLen, h))++opt::nz;
    }
}

void nthashBF(const char *geneName, const char *readName) {

#ifdef _OPENMP
    omp_set_num_threads(opt::threads);
#endif
    std::cerr<<"#threads="<<opt::threads << "\n";
    for(opt::method=0; opt::method<3; opt::method++) {
        std::cerr<<"method="<<itm[opt::method]<<" ";
        for (unsigned k=50; k<=opt::squery; k+=100) {
            opt::kmerLen = k;
            //init_kmod(opt::kmerLen);
            std::cerr<<"kmerl="<<opt::kmerLen<<"\n";
            for (unsigned i=1; i<6; i+=2) {
                opt::nhash = i;
                std::cerr<<"nhash="<<opt::nhash<<" ";
                

                #ifdef _OPENMP
                double sTime = omp_get_wtime();
                #else
                clock_t start = clock();
                #endif

                
                BloomFilter myFilter(opt::ibits*opt::ngene*opt::sgene , opt::nhash, opt::kmerLen);
                loadBf(myFilter, geneName);
                cerr << "|popBF|=" << myFilter.getPop() << " ";

                #ifdef _OPENMP
                std::cerr << "load_time=" <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
                #else
                std::cerr << "load_time=" <<setprecision(4) << (double)(clock() - start)/CLOCKS_PER_SEC << "\n";
                #endif
                
                #ifdef _OPENMP
                sTime = omp_get_wtime();
                #else
                start = clock();
                #endif
                
                queryBf(myFilter, readName);
                
                #ifdef _OPENMP
                std::cerr << "query_time=" <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
                #else
                std::cerr << "query_time=" <<setprecision(4) << (double)(clock() - start)/CLOCKS_PER_SEC << "\n";
                #endif
            }
        }
        cerr << "\n";
    }
}

void nthashRT(const char *readName) {
    getFtype(readName);
    cerr << "CPU time (sec) for hash algorithms for ";
    cerr << "kmer="<<opt::kmerLen<< "\n";
    cerr << "nhash="<<opt::nhash<< "\n";
    for(unsigned method=0; method<2; method++)
        cerr << itm[method] << "\t";
    cerr << "\n";
    if(opt::nhash>1) {
        for(unsigned method=0; method<2; method++) {
            opt::nz=0;
            ifstream uFile(readName);
            string line;
            clock_t sTime = clock();
            while(getSeq(uFile,line)) {
                if(itm[method]=="city")
                    hashSeqcM(line);
                else if(itm[method]=="murmur")
                    hashSeqmM(line);
                else if(itm[method]=="nthash")
                    hashSeqrM(line);
                else if(itm[method]=="ntbase")
                    hashSeqbM(line);
                else if(itm[method]=="xxhash")
                    hashSeqxM(line);
            }
            cerr << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\t";
            uFile.close();
        }
        cerr << "\n";
    }
    else {
        for(unsigned method=0; method<2; method++) {
            opt::nz=0;
            ifstream uFile(readName);
            string line;
            clock_t sTime = clock();
            while(getSeq(uFile, line)) {
                if(itm[method]=="city")
                    hashSeqc(line);
                else if(itm[method]=="murmur")
                    hashSeqm(line);
                else if(itm[method]=="nthash")
                    hashSeqr(line);
                else if(itm[method]=="ntbase")
                    hashSeqb(line);
                else if(itm[method]=="xxhash")
                    hashSeqx(line);
            }
            cerr << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\t";
            uFile.close();
        }
        cerr << "\n";
    }
}

int main(int argc, char** argv) {

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 'j':
            arg >> opt::threads;
            break;
        case 'b':
            arg >> opt::ibits;
            break;
        case 'q':
            arg >> opt::nquery;
            break;
        case 'l':
            arg >> opt::squery;
            break;
        case 't':
            arg >> opt::ngene;
            break;
        case 'g':
            arg >> opt::sgene;
            break;
        case 'h':
            arg >> opt::nhash;
            break;
        case 'k':
            arg >> opt::kmerLen;
            //init_kmod(opt::kmerLen);
            break;
        case 'm':
            arg >> opt::maxitr;
            break;
        case 'a':
            arg >> opt::method;
            break;
        case 'i':
            opt::inpFlag=true;
            break;
        case 'u':
            opt::uniformity=true;
            break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    if (argc - optind != 1 && argc - optind != 2) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }

    if (die) {
        std::cerr << "Try `" << PROGRAM
                  << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

    const char *readName(argv[argc-1]);

    if (opt::uniformity) {
        std::cerr<<"bit/i="<<opt::ibits<<"\n";
        std::cerr<<"nquery="<<opt::nquery<<"\n";
        std::cerr<<"squery="<<opt::squery<<"\n";
        std::cerr<<"ngene="<<opt::ngene<<"\n";
        std::cerr<<"sgene="<<opt::sgene<<"\n\n";

        if(opt::inpFlag) {
            makeGene(opt::ngene, opt::sgene, opt::nquery, opt::squery);
            makeRead(opt::nquery, opt::squery);
            nthashBF("genes.fa", "reads.fa");
        }
        else {
            const char *geneName(argv[argc-2]);
            nthashBF(geneName, readName);
        }
    }
    else {
        if(opt::inpFlag) {
            makeRead(opt::nquery, opt::squery);
            nthashRT("reads.fa");
        }
        else
            nthashRT(readName);
    }

    return 0;
}
