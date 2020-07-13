/**
 * Resolve unitig repeats using a sliding window and
 * and short read information.
 * Written by Vladimir Nikolic <vnikolic@bcgsc.ca>
 */

#include "RAlgorithmsShort.h"

#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#define PROGRAM "rresolver-short"

static const char VERSION_MESSAGE[] =
PROGRAM " (" PACKAGE_NAME ") " VERSION "\n"
"Written by Vladimir Nikolic.\n"
"\n"
"Copyright 2020 Canada's Michael Smith Genome Sciences Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... <contigs> <graph> [<reads1> <reads2> ...]\n"
"Resolve unitig repeats using a sliding window and\n"
"and short read information.\n"
"\n"
" Arguments:\n"
"\n"
"  <contigs>  contigs in FASTA format\n"
"  <graph>    contig adjacency graph\n"
"  <reads>    reads in FASTA format\n"
"\n"
" Options:\n"
"\n"
"  -b, --bloom-size=N          read Bloom filter size.\n"
"                              Unit suffixes 'K' (kilobytes), 'M' (megabytes),\n"
"                              or 'G' (gigabytes) may be used. [required]\n"
"  -g, --graph=FILE            write the contig adjacency graph to FILE. [required]\n"
"  -c, --contigs=FILE          write the contigs to FILE. [required]\n"
"  -j, --threads=N             use N parallel threads [1]\n"
"  -k, --kmer=N                assembly k-mer size\n"
"  -h, --hist=PREFIX           write the algorithm histograms with the given prefix.\n"
"                              Histograms are omitted if no prefix is given."
"  -t, --threshold=N           set path support threshold to N. [4]"
"  -m, --min-tests=N           set minimum number of sliding window moves to N. [20]"
"  -r, --branching=N           set maximum number of branching paths to N. [75]"
"  -e, --error-correction      enable correction of a 1bp error in kmers. [false]"
"  -S, --supported=FILE        write supported paths to FILE.\n"
"  -U, --unsupported=FILE      write unsupported paths to FILE.\n"
"                              Used for path sequence quality check.\n"
"      --adj                   output the graph in ADJ format [default]\n"
"      --asqg                  output the graph in ASQG format\n"
"      --dot                   output the graph in GraphViz format\n"
"      --gfa                   output the graph in GFA1 format\n"
"      --gfa1                  output the graph in GFA1 format\n"
"      --gfa2                  output the graph in GFA2 format\n"
"      --gv                    output the graph in GraphViz format\n"
"      --sam                   output the graph in SAM format\n"
"  -v, --verbose               display verbose output\n"
"      --help                  display this help and exit\n"
"      --version               output version information and exit\n"
"\n"
"Report bugs to <" PACKAGE_BUGREPORT ">.\n";

namespace opt {

  /** Read Bloom filter size in bytes. */
  size_t bloomSize = 0;

  /** The number of parallel threads */
  int threads = 1;

  /** Prefix for the histogram files */
  std::string histPrefix;

  /** Name of the file to write resulting graph to */
  std::string outputGraphPath;

  /** Name of the file to write resulting graph to */
  std::string outputContigsPath;

  /** Number of kmers required to be found for a path to be supported */
  int threshold = 4;

  /** Minimum number of sliding window moves */
  int min_tests = 20;

  /** Maximum number of branching paths */
  int branching = 75;

  /** Flag indicating whether error correction is enabled */
  int error_correction = 0;

  /** Name of the file to write supported paths to */
  std::string outputSupportedPathsPath;

  /** Name of the file to write unsupported paths to */
  std::string outputUnsupportedPathsPath;

  unsigned k = 0;  // used by ContigProperties

  int format;  // used by ContigProperties

}

static const char shortopts[] = "b:j:g:c:k:h:t:m:r:eS:U:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
  { "bloom-size",        required_argument, NULL, 'b' },
  { "threads",           required_argument, NULL, 'j' },
  { "graph",             required_argument, NULL, 'g' },
  { "contigs",           required_argument, NULL, 'c' },
  { "kmer",              required_argument, NULL, 'k' },
  { "hist",              required_argument, NULL, 'h' },
  { "threshold",         required_argument, NULL, 't' },
  { "min-tests",         required_argument, NULL, 'm' },
  { "branching",         required_argument, NULL, 'r' },
  { "error-correction",  no_argument,       &opt::error_correction, 1},
  { "supported",         required_argument, NULL, 'S' },
  { "unsupported",       required_argument, NULL, 'U' },
  { "adj",               no_argument,       &opt::format, ADJ  },
  { "asqg",              no_argument,       &opt::format, ASQG },
  { "dot",               no_argument,       &opt::format, DOT  },
  { "gfa",               no_argument,       &opt::format, GFA1 },
  { "gfa1",              no_argument,       &opt::format, GFA1 },
  { "gfa2",              no_argument,       &opt::format, GFA2 },
  { "gv",                no_argument,       &opt::format, DOT  },
  { "sam",               no_argument,       &opt::format, SAM  },
  { "verbose",           no_argument,       NULL, 'v' },
  { "help",              no_argument,       NULL, OPT_HELP },
  { "version",           no_argument,       NULL, OPT_VERSION },
  { NULL, 0, NULL, 0 }
};

void check_options(int argc, bool& die) {
  if (opt::bloomSize == 0) {
    std::cerr << PROGRAM ": missing or invalid value for mandatory option `-b'" << std::endl;
    die = true;
  }

  if (argc - optind < 3) {
    std::cerr << PROGRAM ": missing input file arguments" << std::endl;
    die = true;
  }

  if (opt::k <= 0) {
    std::cerr << PROGRAM ": missing or invalid value for mandatory option `-k'" << std::endl;
    die = true;
  }

  if (opt::outputGraphPath.empty()) {
    std::cerr << PROGRAM ": missing or invalid value for mandatory option `-g`" << std::endl;
    die = true;
  }

  if (opt::outputContigsPath.empty()) {
    std::cerr << PROGRAM ": missing or invalid value for mandatory option `-c`" << std::endl;
    die = true;
  }

  if (opt::threads <= 0) {
    std::cerr << PROGRAM ": invalid number of threads `-j`" << std::endl;
    die = true;
  }
}

void writePaths(std::ostream& stream, const ImaginaryContigPaths& paths)
{
  int counter = 0;
  for (const auto &path : paths) {
    stream << '>' << counter++ << '\n' << getPathSequence(path) << '\n';
  }
}

void writeResults(const ImaginaryContigPaths &supportedPaths,
                  const ImaginaryContigPaths &unsupportedPaths,
                  const std::string &commandLine)
{
  if (!opt::outputContigsPath.empty()) {
    storeContigs(opt::outputContigsPath);
	}

  if (!opt::outputGraphPath.empty()) {
    storeContigGraph(opt::outputGraphPath, PROGRAM, commandLine);
	}

  if (!opt::outputSupportedPathsPath.empty()) {
    if (opt::verbose) {
      std::cerr << "Writing supported paths to `" << opt::outputSupportedPathsPath << "'..." << std::endl;
    }
    std::ofstream outputSupportedPaths(opt::outputSupportedPathsPath);
    writePaths(outputSupportedPaths, supportedPaths);
    assert(outputSupportedPaths.good());
    outputSupportedPaths.close();
    if (opt::verbose) {
      std::cerr << "Supported paths written." << std::endl;
    }
  }

  if (!opt::outputUnsupportedPathsPath.empty()) {
    if (opt::verbose) {
      std::cerr << "Writing unsupported paths to `" << opt::outputUnsupportedPathsPath << "'..." << std::endl;
    }
    std::ofstream outputUnsupportedPaths(opt::outputUnsupportedPathsPath);
    writePaths(outputUnsupportedPaths, unsupportedPaths);
    assert(outputUnsupportedPaths.good());
    outputUnsupportedPaths.close();
    if (opt::verbose) {
      std::cerr << "Unsupported paths written." << std::endl;
    }
  }
}

int main(int argc, char** argv)
{
	std::string commandLine;
	{
		std::ostringstream ss;
		char** last = argv + argc - 1;
		std::copy(argv, last, std::ostream_iterator<const char *>(ss, " "));
		ss << *last;
		commandLine = ss.str();
	}

  bool die = false; //scary, and spooky
  for (int c; (c = getopt_long(argc, argv,
          shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
      case '?':
        die = true; break;
      case 'b':
        opt::bloomSize = SIToBytes(arg); break;
      case 'j':
        arg >> opt::threads; break;
      case 'k':
        arg >> opt::k; break;
      case 'h':
        arg >> opt::histPrefix; break;
      case 'g':
        arg >> opt::outputGraphPath; break;
      case 'c':
        arg >> opt::outputContigsPath; break;
      case 't':
        arg >> opt::threshold; break;
      case 'm':
        arg >> opt::min_tests; break;
      case 'r':
        arg >> opt::branching; break;
      case 'S':
        arg >> opt::outputSupportedPathsPath; break;
      case 'U':
        arg >> opt::outputUnsupportedPathsPath; break;
      case 'v':
        ++opt::verbose; break;
      case OPT_HELP:
        std::cout << USAGE_MESSAGE; exit(EXIT_SUCCESS);
      case OPT_VERSION:
        std::cout << VERSION_MESSAGE; exit(EXIT_SUCCESS);
    }

    if (optarg != NULL && (!arg.eof() || arg.fail())) {
      std::cerr << PROGRAM ": invalid option: `-"
        << char(c) << optarg << "'" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  check_options(argc, die);

  if (die) {
    std::cerr << "Try `" << PROGRAM
      << " --help' for more information." << std::endl;
    exit(EXIT_FAILURE);
  }

#if _OPENMP
  if (opt::threads > 0)
    omp_set_num_threads(opt::threads);
#endif

  std::string contigsPath(argv[optind++]);
  std::string contigGraphPath(argv[optind++]);

  loadContigGraph(contigGraphPath);
  loadContigs(contigsPath);

  std::vector<std::string> readFilepaths;
  for (int i = optind; i < argc; i++) {
    readFilepaths.push_back(argv[i]);
  }

  ImaginaryContigPaths supportedPaths, unsupportedPaths;
  resolveShort(readFilepaths, supportedPaths, unsupportedPaths);
  
  if (opt::verbose) {
    std::cerr << "Stats after resolution:" << std::endl;
    printGraphStats(std::cerr, g_contigGraph);
  }

  writeResults(supportedPaths, unsupportedPaths, commandLine);

  return EXIT_SUCCESS;
}
