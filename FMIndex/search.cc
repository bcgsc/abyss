/*
 *  Copyright (c) 2010 Yasuo Tabei
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
*/

#include "FMIndex.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int mode = 0; // exact:0 hamming:1 edit:2
int dist = 0; // distance threshold

static void usage();
static void parse_parameters (int argc, char **argv);

int main(int argc, char **argv)
{
	parse_parameters(argc, argv);
	FMIndex f;
	ifstream is(argv[argc-2]);
	f.load(is);
	vector<vector<uint8_t> > qs;
	f.readQstring(argv[argc-1], qs);
	double sTime = clock();
	if (mode == 0) {
		cerr << "exact search mode:" << endl;
		pair<uint64_t, uint64_t> res;
		for (size_t i = 0; i < qs.size(); i++) {
			vector<uint8_t> &s = qs[i];
			f.search(s, res);
			cout << "query id: " << i << endl;
			if (res.first != (size_t)-1 || res.first != (size_t)-1) {
				for (size_t it = res.first; it <= res.second; it++)
					cout << f.locate(it) << " ";
				cout << endl;
			}
		}
	}
	else if (mode == 1) {
		cerr << "Hamming distance mode:" << endl;
		for (size_t i = 0; i < qs.size(); i++) {
			vector<uint8_t> &s = qs[i];
			vector<pair<uint64_t, uint64_t> > res;
			f.searchHamming(s, dist, res);
			cout << "query id:" << i << endl;
			for (size_t i = 0; i < res.size(); i++) {
				for (size_t it = res[i].first; it <= res[i].second; it++)
					cout << f.locate(it) << " ";
			}
			cout << endl;
		}
	}
	else if (mode == 2) {
		cerr << "Edit distance mode:" << endl;
		for (size_t i = 0; i < qs.size(); i++) {
			vector<uint8_t> &s = qs[i];
			vector<pair<uint64_t, uint64_t> > res;
			f.searchEdit(s, dist, res);
			cout << "query id:" << i << endl;
			for (size_t i = 0; i < res.size(); i++) {
				for (size_t it = res[i].first; it <= res[i].second; it++)
					cout << f.locate(it) << " ";
			}
			cout << endl;
		}
	}

	double eTime = clock();
	fprintf(stderr, "cpu time: %f\n", (eTime - sTime)/CLOCKS_PER_SEC);

	return 0;
}

static void usage()
{
	std::cerr << std::endl
		<< "Usage: fmsearch [OPTION]... INDEXFILE INFILE" << std::endl << std::endl
		<< "       where [OPTION]...  is a list of zero or more optional arguments" << std::endl
		<< "             INDEXFILE    is the name of index-file" << std::endl
		<< "             INFILE       is the name of an input file" << std::endl
		<< "             -mode [0|1|2]" << std::endl
		<< "             1:exact search 2:hamming distance search 3:edit distance search" << std::endl
		<< "             (default: " << mode << ")" << std::endl
		<< "             -dist [maximum hamming or edit distance]" << std::endl
		<< "             (default: " << dist << ")" << std::endl
		<< std::endl;
	exit(0);
}

static void parse_parameters (int argc, char **argv)
{
	if (argc == 1) usage();
	int argno;
	for (argno = 1; argno < argc-1; argno++){
		if (argv[argno][0] == '-'){
			if (!strcmp (argv[argno], "-mode")) {
				if (argno == argc - 1) std::cerr << "Must specify miximum itemset size after -hamdist" << std::endl;
				mode = atoi(argv[++argno]);
			}
			else if (!strcmp (argv[argno], "-dist")) {
				if (argno == argc - 1) std::cerr << "Must specify miximum itemset size after -dist" << std::endl;
				dist = atoi(argv[++argno]);
			}
			else {
				usage();
			}
		} else {
			break;
		}
	}
}
