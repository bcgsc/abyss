#ifndef SEQGEN_H_
#define SEQGEN_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

const int alphNum = 4;
const char iTb[alphNum] = {'A', 'C', 'G', 'T'};

size_t sequenceLen;

void makeRead(const size_t rNum, const size_t rLen) {
    srand(time(NULL));
    size_t bTc[alphNum] = {0, 0, 0, 0};
    string myrSeq;
    myrSeq.resize(rLen);
    ofstream rFile("reads.fa");

    for (size_t j=0; j < rNum; j++) {
        for (size_t i=0; i< rLen; i++) {
            int ranInd = rand() % alphNum;
            myrSeq[i] = iTb[ranInd];
            ++bTc[ranInd];
        }
        rFile << ">f" << j << "\n" << myrSeq << "\n";
    }

    //for (size_t i=0; i< alphNum; i++)
    //    cerr << bTc[i] << "\n";
    rFile.close();
}


void makeGenome(const size_t seqLen) {
    srand(time(NULL));
    size_t bTc[alphNum] = {0, 0, 0, 0};
    string mygSeq;
    mygSeq.resize(seqLen);
    for (size_t i=0; i< seqLen; i++) {
        int ranInd = rand() % alphNum;
        mygSeq[i] = iTb[ranInd];
        ++bTc[ranInd];
    }
    ofstream gFile("genome.fa");
    gFile << ">1\n" << mygSeq << "\n";
    //for (int i=0; i< alphNum; i++)
    //    cerr << bTc[i] << "\n";
    gFile.close();
}

void makeGene(const size_t gNum, const size_t gLen) {
    srand(time(NULL));
    size_t bTc[alphNum] = {0, 0, 0, 0};
    string mygSeq;
    mygSeq.resize(gLen);
    ofstream gFile("genes.fa");
    for (size_t j=0; j < gNum; j++) {
        for (size_t i=0; i< gLen; i++) {
            int ranInd = rand() % alphNum;
            mygSeq[i] = iTb[ranInd];
            ++bTc[ranInd];
        }
    }
    gFile.close();
}

void makeGene(const size_t gNum, const size_t gLen, const size_t rNum, const size_t rLen) {
    srand(time(NULL));
    size_t bTc[alphNum] = {0, 0, 0, 0};
    string mygSeq;
    mygSeq.resize(gLen);
    ofstream gFile("genes.fa");
    ofstream rFile("treads.fa");
    for (size_t j=0; j < gNum; j++) {
        for (size_t i=0; i< gLen; i++) {
            int ranInd = rand() % alphNum;
            mygSeq[i] = iTb[ranInd];
            ++bTc[ranInd];
        }
        gFile << ">g" << j << "\n" << mygSeq << "\n";
        for(size_t k=0; k<rNum/gNum; k++) {
            int ranR = rand() % (gLen-rLen+1);
            rFile << ">t" << k << "\n" << mygSeq.substr(ranR, rLen) << "\n";
        }

    }
    gFile.close();
    rFile.close();
}


int testseqgen() {
    //makeGenome(3000000000+32-1);
    makeRead(1000000000, 32);
    makeGene(50, 2000000+32-1);
    //makeGene(6, 50+32-1);
    return 0;
}

#endif