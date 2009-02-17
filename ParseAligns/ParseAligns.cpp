#include "Aligner.h"
#include "PairUtils.h"
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

using namespace std;

// TYPEDEFS
typedef std::map<std::string, AlignmentVector> ReadAlignMap;


// FUNCTIONS
bool checkUniqueAlignments(int kmer, const AlignmentVector& alignVec);

std::string makePairID(std::string refID);

/**
 * Return the size of the fragment demarcated by the specified
 * alignments.
 * @return INT_MIN if the pair is invalid
 */
static int fragmentSize(const Alignment& a0, const Alignment& a1)
{
	assert(a0.contig == a1.contig);
	if (a0.isRC == a1.isRC)
		return INT_MIN;
	const Alignment& a = a0.isRC ? a1 : a0;
	const Alignment& b = a0.isRC ? a0 : a1;
	return b.targetAtQueryEnd() - a.targetAtQueryStart();
}

// Read in the alignments file into the table
static void readAlignments(istream& in, ReadAlignMap* pout)
{
	ReadAlignMap& out = *pout;
	string line;
	for (string line; getline(in, line);) {
		stringstream s(line);
		string readID;
		s >> readID;
		AlignmentVector& alignments = out[readID];
		assert(alignments.empty());
		for (Alignment ali; s >> ali;)
			alignments.push_back(ali);
	}
	assert(in.eof());
}

static void assert_open(ifstream& f, const string& p)
{
	if (!f.is_open()) {
		cerr << p << ": " << strerror(errno) << endl;
		exit(EXIT_FAILURE);
	}
}

static void readAlignmentsFile(string path, ReadAlignMap* pout)
{
	cout << "Reading `" << path << "'..." << endl;
	ifstream fin(path.c_str());
	assert_open(fin, path);
	readAlignments(fin, pout);
	fin.close();
}

int main(int argc, char* const* argv)
{

	if(argc < 2)
	{
		std::cout << "Usage: <kmer size> <list of files to parse>\n";
		exit(1);
	}

	int kmer = atoi(argv[1]);
	assert(kmer > 0);
	int optind = 2;

	std::ofstream pairedAlignFile("PairAligns.txt");
	std::ofstream distanceList("DistanceList.txt");

	ReadAlignMap alignTable;

	if (optind < argc) {
		for_each(argv + optind, argv + argc,
				bind2nd(ptr_fun(readAlignmentsFile), &alignTable));
	} else {
		cout << "Reading from standard input..." << endl;
		readAlignments(cin, &alignTable);
	}

	std::cout << "Align table has " << alignTable.size() << " entries\n";

	int numDifferent = 0;
	int numSame = 0;
	int numInvalid = 0;
	int numMissed = 0;
	int numMulti = 0;
	int numNonSingle = 0;
	
	// parse the alignment table
	for(ReadAlignMap::iterator iter = alignTable.begin(); iter != alignTable.end(); ++iter)
	{
		std::string currID = iter->first;
		std::string pairID = makePairID(currID);
		
		// Find the pair align
		ReadAlignMap::iterator pairIter = alignTable.find(pairID);
		if(pairIter != alignTable.end())
		{
			// Both reads must align to a unique location.
			// The reads are allowed to span more than one contig, but
			// at least one of the two reads must span no more than
			// two contigs.
			bool isRefUnique = checkUniqueAlignments(kmer, iter->second);
			bool isPairUnique = checkUniqueAlignments(kmer, pairIter->second);
			const unsigned MAX_SPAN = 2;
			if ((iter->second.size() <= MAX_SPAN
						|| pairIter->second.size() <= MAX_SPAN)
					&& isRefUnique && isPairUnique) {
				// Iterate over the vectors, outputting the aligments
				for(AlignmentVector::iterator refAlignIter = iter->second.begin(); refAlignIter != iter->second.end(); ++refAlignIter)
				{
					for(AlignmentVector::iterator pairAlignIter = pairIter->second.begin(); pairAlignIter != pairIter->second.end(); ++pairAlignIter)
					{
						// Are they on the same contig and the ONLY alignments?
						if(refAlignIter->contig == pairAlignIter->contig)
						{
							if((iter->second.size() == 1 && pairIter->second.size() == 1))
							{
								int size = fragmentSize(
										*refAlignIter,
										*pairAlignIter);
								if (size > INT_MIN) {
									distanceList << size << "\n";
									numSame++;
								} else
									numInvalid++;
							}
						}
						else
						{
							// Print the alignment and the swapped alignment
							pairedAlignFile << currID << " " << *refAlignIter << " " << pairID << " " << *pairAlignIter << "\n";
							pairedAlignFile << pairID << " " << *pairAlignIter << " " << currID << " " << *refAlignIter << "\n";
							numDifferent++;
						}							
					}
				}
				
				
			}
			else
			{
				if(!isRefUnique || !isPairUnique)
				{
					numMulti++;
				}
				else
				{
					numNonSingle++;
				}
			}

			
			// Erase the pair as its not needed (faster to mark it as invalid?)
			alignTable.erase(pairIter);
		}
		else
		{
			numMissed++;
		}
	}
	
	printf("Num unmatched: %d Num same contig: %d Invalid: %d "
			"Num diff contig: %d Num multi: %d "
			"Num not singular: %d\n",
			numMissed, numSame, numInvalid, numDifferent, numMulti,
			numNonSingle);
	
	pairedAlignFile.close();
	distanceList.close();
}

bool checkUniqueAlignments(int kmer, const AlignmentVector& alignVec)
{
	// Ensure that no read start position hit to more than 1 contig
	assert(!alignVec.empty());
	
	const int num_starts = alignVec.front().read_length - kmer + 1;
	int* coverage = new int[num_starts];
	
	for(int i = 0; i < num_starts; ++i)
	{
		coverage[i] = 0;
	}
	
	for(AlignmentVector::const_iterator iter = alignVec.begin(); iter != alignVec.end(); ++iter)
	{
		//std::cout << "Align: " << *iter << "\n";
		int length = iter->align_length;
		int start = iter->read_start_pos;
		for(int i = 0; i < (length - kmer + 1); ++i)
		{
			int start_idx = start + i;
			assert(start_idx >= 0 && start_idx < num_starts);
			coverage[start_idx]++;
		}
	}
	
	bool unique = true;
	//printf("Coverage: \n");
	for(int i = 0; i < num_starts; ++i)
	{
		//printf("%d", coverage[i]);
		if(coverage[i] > 1)
		{
			unique = false;
		}
	}
	//printf("\n");
	delete [] coverage;
	return unique;
}

// Change the id into the id of its pair
std::string makePairID(std::string refID)
{
	std::string pairID = refID;
	// Change the last character
	ssize_t lastIdx = pairID.size() - 1;
	assert(lastIdx > 0);
	char c = refID[lastIdx];
	bool matched = true;
	switch (c) {
		case '1': c = '2'; break;
		case '2': c = '1'; break;
		case 'A': c = 'B'; break;
		case 'B': c = 'A'; break;
		case 'F': c = 'R'; break;
		case 'R': c = 'F'; break;
		default: matched = false; break;
	}
	if (matched) {
		pairID[lastIdx] = c;
	} else {
		static const char *suffix0 = "forward";
		static const char *suffix1 = "reverse";
		ssize_t i = refID.length() - strlen(suffix0);
		assert(i > 0);
		string suffix = refID.substr(i);
		if (suffix == suffix0) {
			pairID.replace(i, string::npos, suffix1);
		} else if (suffix == suffix1) {
			pairID.replace(i, string::npos, suffix0);
		} else {
			cerr << "error: read ID `" << refID
				<< "' must end in one of\n"
				"\t/1 and /2 or"
				" _A and _B or"
				" _F and _R or"
				" _forward and _reverse\n";
			exit(EXIT_FAILURE);
		}
	}
	return pairID;
}
