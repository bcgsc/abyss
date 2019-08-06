#include "../BloomFilterUtil.h"
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>

using namespace std;

int
getFileSize(string file)
{
	ifstream is(file.c_str(), ifstream::ate | ifstream::binary);
	return is.tellg();
}

void
containSeq(BloomFilter& bloom, const string& seq, unsigned hashNum, unsigned kmerSize)
{
	int total = 0;
	int hits = 0;
	for (int i = 0; i < seq.length() - kmerSize; i++) {
		total++;
		string kmer = seq.substr(i, kmerSize);
		if (bloom.contains(kmer.c_str()))
			hits++;
		else
			cout << "Missing: " << kmer << endl;
	}
	cerr << "Found " << hits << " out of " << total << endl;
}

void
contigsToBloom(const string& file, BloomFilter& bloom, unsigned hashNum, unsigned kmerSize)
{
	ifstream infile(file.c_str());
	clock_t t;

	if (!infile) {
		cerr << "Can't open file " << file << std::endl;
		exit(EXIT_FAILURE);
	}

	cout << "Contigs processed k=" << kmerSize << ":\n";

	string line;
	string seq = "";
	string prevHead = "";
	int cttig = 0;

	while (getline(infile, line)) {
		if (line.compare(0, 1, ">") == 0) {
			if (line.compare(prevHead) != 0 && seq.length() != 0 && prevHead.length() != 0) {
				cttig++;
				cout << cttig << endl;
				containSeq(bloom, seq, hashNum, kmerSize);
			}
			seq = "";
			prevHead = line;
		} else {
			seq += line;
		}
	}
	cttig++;
	cout << cttig << endl;
	containSeq(bloom, seq, hashNum, kmerSize);
	infile.close();
}

int
main(int argc, const char* argv[])
{
	string myFile;
	double fpr = 0.0001;
	unsigned kmerSize = 15;

	if (argc < 3) {
		cerr << "Usage: " << endl;
		cerr << "-f  sequences to scaffold (Multi-FASTA format, required)" << endl;
		cerr << "-k  k-mer value (default -k " << kmerSize << ", optional)" << endl;
		cerr << "-p  Bloom filter false positive rate (default -p " << fpr
		     << " optional - increase to prevent memory allocation errors)" << endl;
		return 1;
	}

	for (int i = 1; i < argc; i++) {
		string arg = argv[i];
		if (arg == "-f") {
			if (i + 1 < argc) {
				myFile = argv[++i];
			} else
				cerr << "Must have an argument after the -f option" << endl;
		} else if (arg == "-k") {
			if (i + 1 < argc) {
				kmerSize = atoi(((string)argv[++i]).c_str());
			}
		} else if (arg == "-p") {
			if (i + 1 < argc) {
				fpr = atof(((string)argv[++i]).c_str());
			}
		}
	}

	if (myFile.size() == 0) {
		cerr << "Must include the -f option" << endl;
		return 1;
	}

	int bfelements = getFileSize(myFile);
	size_t size = ceil((-1 * bfelements * log(fpr)) / (log(2) * log(2)));
	size += 64 - (size % 64);
	unsigned hashNum = floor((size / bfelements) * log(2));
	cout << "***** Bloom filter specs: \nelements = " << bfelements << "\nFPR = " << fpr
	     << "\nsize (bits) = " << size << "\nhash functions = " << hashNum << endl;

	BloomFilter bloom("BloomFilter.bf");

	cout << "Shredding supplied sequence file (-f " << myFile << ") into " << kmerSize
	     << "-mers..\n";

	contigsToBloom(myFile, bloom, hashNum, kmerSize);

	cout << "Done!" << endl;

	return 0;
}
