#include <stdio.h>
#include <iostream>
#include <stdlib.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

void printBin(uint64_t n) {
    int count=0,count1=0;
    while (++count<=64) {
        if (n & ((uint64_t)1<<63)) {
            printf("1 ");
            ++count1;
        }
        else
            printf("0 ");

        n <<= 1;
    }
    printf(" count1=%d\n",count1);
}

void balanceSeed() {

    srand(time(NULL));
    int r1[4]= {32,32,32,32};
    const int seedNum=4;
    int hashSeed[seedNum][64];
    for (int i=0; i<seedNum; i++)
        for (int j=0; j<64; j++)
            hashSeed[i][j]=0;
    for (int j=0; j<64; j++) {
        int ranVec[seedNum];
        for (int i=0; i<seedNum; i++)
            ranVec[i]=i;
        for (int i=0; i<seedNum/2; i++) {
            int ranInd = rand() % (seedNum-i);
            int tmp = ranVec[seedNum-i-1];
            ranVec[seedNum-i-1] = ranVec[ranInd];
            ranVec[ranInd] = tmp;
        }
        if(r1[ranVec[0]]>0 && r1[ranVec[1]]>0) {
            hashSeed[ranVec[0]][j] = 1;
            r1[ranVec[0]]--;
            hashSeed[ranVec[1]][j] = 1;
            r1[ranVec[1]]--;
        }
        else if(r1[ranVec[0]]<0 && r1[ranVec[1]]>0) {
            hashSeed[ranVec[1]][j] = 1;
            r1[ranVec[1]]--;
            if(r1[ranVec[2]]>0) {
                hashSeed[ranVec[2]][j] = 1;
                r1[ranVec[2]]--;
            }
            else
                hashSeed[ranVec[3]][j] = 1;
            r1[ranVec[3]]--;
        }
        else if(r1[ranVec[0]]>0 && r1[ranVec[1]]<0) {
            hashSeed[ranVec[0]][j] = 1;
            r1[ranVec[0]]--;
            if(r1[ranVec[2]]>0) {
                hashSeed[ranVec[2]][j] = 1;
                r1[ranVec[2]]--;
            }
            else
                hashSeed[ranVec[3]][j] = 1;
            r1[ranVec[3]]--;
        }
        else {
            hashSeed[ranVec[2]][j] = 1;
            r1[ranVec[2]]--;
            hashSeed[ranVec[3]][j] = 1;
            r1[ranVec[3]]--;
        }

    }

    for (int i=0; i<seedNum; i++) {
        for (int j=0; j<64; j++)
            std::cout << hashSeed[i][j] << " ";
        std::cout << "\n";
    }

    long long unsigned int allSeed[seedNum];
    for (int i=0; i<seedNum; i++) {
        uint64_t hSeed=0;
        for (int j=0; j<64; j++) {
            if(hashSeed[i][j]==1)
                hSeed = (uint64_t) (hSeed << 1 | 1);
            else
                hSeed = (uint64_t) (hSeed << 1 | 0);
        }
        allSeed[i]=hSeed;
        printf("%" PRIu64 "\n", hSeed);
        //printf("0x%llx\n", hSeed);
        printBin(hSeed);

    }
    //for (int i=0; i<seedNum; i++)
    printf("static const uint64_t seedA = 0x%llx;\n", allSeed[0]);
    printf("static const uint64_t seedC = 0x%llx;\n", allSeed[1]);
    printf("static const uint64_t seedG = 0x%llx;\n", allSeed[2]);
    printf("static const uint64_t seedT = 0x%llx;\n", allSeed[3]);
    exit(0);
}

void getSeeds() {
    //balanceSeed();
    srand(time(NULL));
    const int seedNum=4;
    int hashSeed[seedNum][64];
    for (int i=0; i<seedNum; i++)
        for (int j=0; j<64; j++)
            hashSeed[i][j]=0;


    for (int j=0; j<64; j++) {
        int ranVec[seedNum];
        for (int i=0; i<seedNum; i++)
            ranVec[i]=i;
        for (int i=0; i<seedNum/2; i++) {
            int ranInd = rand() % (seedNum-i);
            int tmp = ranVec[seedNum-i-1];
            ranVec[seedNum-i-1] = ranVec[ranInd];
            ranVec[ranInd] = tmp;
        }
        for (int i=0; i<seedNum/2; i++)
            hashSeed[ranVec[i]][j] = 1;
    }

    for (int i=0; i<seedNum; i++) {
        for (int j=0; j<64; j++)
            std::cout << hashSeed[i][j] << " ";
        std::cout << "\n";
    }

    long long unsigned int allSeed[seedNum];

    for (int i=0; i<seedNum; i++) {
        uint64_t hSeed=0;
        for (int j=0; j<64; j++) {
            if(hashSeed[i][j]==1)
                hSeed = (uint64_t) (hSeed << 1 | 1);
            else
                hSeed = (uint64_t) (hSeed << 1 | 0);
        }
        allSeed[i]=hSeed;
        printf("%" PRIu64 "\n", hSeed);
        printBin(hSeed);

    }

    printf("static const uint64_t seedA = 0x%llx;\n", allSeed[0]);
    printf("static const uint64_t seedC = 0x%llx;\n", allSeed[1]);
    printf("static const uint64_t seedG = 0x%llx;\n", allSeed[2]);
    printf("static const uint64_t seedT = 0x%llx;\n", allSeed[3]);
    //return 0;
}
