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
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

using namespace std;

char *fname, *oname;
int percent = 25;

void usage();
void parse_parameters (int argc, char **argv);

int main(int argc, char **argv) {
  parse_parameters(argc, argv);

  FMIndex f;
  f.buildFmIndex(fname, percent);

  ofstream os(oname);
  f.save(os);

  return 0;
}

void usage(){
  std::cerr << std::endl
       << "Usage: fmconstruct [OPTION]... INFILE INDEXFILE" << std::endl << std::endl
       << "       where [OPTION] is a list of zero or more optional arguments" << std::endl
       << "             INFILE       is the name of an input file" << std::endl
       << "             INDEXFILE    is the name of index-file" << std::endl
       << "Additional arguments (input and output files may be specified):"  << std::endl
       << "             -percent [val: 0<= val <= 100]: percentage of sampling suffix array positions" << std::endl
       << "             (default: " <<  percent << ")" << std::endl
       << std::endl;
  exit(0);
}

void parse_parameters (int argc, char **argv){
  if (argc == 1) usage();
  int argno;
  for (argno = 1; argno < argc; argno++){
    if (argv[argno][0] == '-'){
      if (!strcmp (argv[argno], "-percent")) {
	if (argno == argc - 1) std::cerr << "Must specify a float value after -percent" << std::endl;
	percent = atof(argv[++argno]);
	if (percent < 0 || percent > 100) {
	  usage();
	}
      }
      else {
	usage();
      }
    } else {
      break;
    }
  }
  if (argno > argc)
    usage();

  fname = argv[argno];
  oname = argv[argno+1];
}
