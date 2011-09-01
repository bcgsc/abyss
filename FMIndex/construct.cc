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
#include "IOUtil.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

using namespace std;

char *fname, *oname;
static unsigned g_sampleSA = 4;

void parse_parameters (int argc, char **argv);

int main(int argc, char **argv) {
  parse_parameters(argc, argv);

  FMIndex fm;
  cerr << "Reading `" << fname << "'...\n";
  {
	  vector<FMIndex::value_type> s;
	  readFile(fname, s);
	  fm.assign(s.begin(), s.end());
  }
  fm.sampleSA(g_sampleSA);

  ofstream out(oname);
  out << fm;
  assert_good(out, oname);

  return 0;
}

/** Print the usage message. */
static void usage()
{
	cout <<
		"Usage: fmconstruct [OPTION]... FILE INDEX\n"
		"Index FILE and store the index in INDEX.\n"
		"\n"
		"  -sample N  period of sampling the suffix array"
			" [" << g_sampleSA << "]\n";
	exit(EXIT_SUCCESS);
}

void parse_parameters (int argc, char **argv){
  if (argc < 3)
	  usage();
  int argno;
  for (argno = 1; argno < argc; argno++){
    if (argv[argno][0] == '-'){
      if (!strcmp (argv[argno], "-sample")) {
	if (argno == argc - 1) std::cerr << "Must specify an integer value after -sample" << std::endl;
	g_sampleSA = atoi(argv[++argno]);
	if (g_sampleSA == 0)
	  usage();
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
