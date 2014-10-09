#include "AssemblyAlgorithms.h"
#include "Common/InsOrderedMap.h"
#include <string>
#include <vector>

namespace AssemblyAlgorithms
{
static inline void loadSequences(ISequenceCollection* seqCollection,std::string inFile)
{
	Timer timer("LoadSequences " + inFile);

	logger(0) << "Reading `" << inFile << "'...\n";
	size_t count = 0, count_good = 0,
	count_small = 0, count_nonACGT = 0,
	count_reversed = 0;
	int fastaFlags = opt::maskCov ?  FastaReader::NO_FOLD_CASE :FastaReader::FOLD_CASE;
	FastaReader reader(inFile.c_str(), fastaFlags);
	
	for (FastaRecord rec; reader >> rec;) 
	{
		//rec sequence loaded into seq object
		Sequence seq = rec.seq;
		size_t len = seq.length();
		if ((opt::kmerSize * 2) > len) 
		{
			count_small++;
			continue;
		}
		//color-space is handled - unnecessary
		if (opt::rank <= 0
				&& count == 0 && seqCollection->empty()) {
			// Detect colour-space reads.
			bool colourSpace
				= seq.find_first_of("0123") != string::npos;
			seqCollection->setColourSpace(colourSpace);
			if (colourSpace)
				cout << "Colour-space assembly\n";
		}
		if (isalnum(seq[0])) 
		{
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));
		}	
		bool good = seq.find_first_not_of("ACGT0123") == string::npos;		
		bool discarded = true;
		if (opt::ss && rec.id.size() > 2 && rec.id.substr(rec.id.size()-2) == "/1") 
		{
			seq = reverseComplement(seq);
			count_reversed++;
		}
		for (unsigned i = 0; i < len - opt::kmerSize - opt::delta + 1; i++) 
		{
			Sequence kmer;
			if (opt::delta > 0) 
			{
				stringstream sseed;
				sseed << Sequence(seq, i, opt::kmerSize) << Sequence(seq, i + opt::kmerSize + opt::delta, opt::kmerSize);
				kmer = sseed.str();
				if (sseed.str().size() != (opt::kmerSize*2)) 
				{
					cout << Sequence(seq, i, opt::kmerSize + opt::delta)
						<< '\n'
						<< Sequence(seq, i, opt::kmerSize) << '\n'
						<< Sequence(seq, i + opt::kmerSize + opt::delta, opt::kmerSize)
						<< '\n';
					exit(EXIT_FAILURE);
				}
			} 
			// if delta=0
			else
				kmer = Sequence(seq, i, opt::kmerSize);
			
			if (good || kmer.find_first_not_of("acgtACGT0123") == string::npos) 
			{
				if (good || kmer.find_first_of("acgt") == string::npos)
				{
					Sequence buf1(kmer,0,opt::kmerSize);
					Sequence buf2(kmer,opt::kmerSize,opt::kmerSize);
					seqCollection->add(KmerPair(buf1,buf2));	
				}
				else 
				{
					transform(kmer.begin(), kmer.end(), kmer.begin(),::toupper);
					Sequence buf1(kmer,0,opt::kmerSize);
					Sequence buf2(kmer,opt::kmerSize,opt::kmerSize);
					seqCollection->add(KmerPair(buf1,buf2), 0);
				}
				discarded = false;	
			}
		}
		
		if (discarded)
			count_nonACGT++;
		else
			count_good++;

		if (++count % 100000 == 0) 
		{
			logger(1) << "Read " << count << " reads. ";
			seqCollection->printLoad();
		}
		seqCollection->pumpNetwork();
	}
	assert(reader.eof());

	logger(1) << "Read " << count << " reads. ";
	seqCollection->printLoad();

	if (count_reversed > 0)
		cerr << "`" << inFile << "': "
			"reversed " << count_reversed << " reads\n";
	if (count_small > 0)
		cerr << "`" << inFile << "': "
			"discarded " << count_small << " reads "
			"shorter than " << opt::kmerSize << " bases\n";
	if (reader.unchaste() > 0)
		cerr << "`" << inFile << "': "
			"discarded " << reader.unchaste() << " unchaste reads\n";
	if (count_nonACGT > 0)
		cerr << "`" << inFile << "': "
			"discarded " << count_nonACGT << " reads "
			"containing non-ACGT characters\n";
	if (count_good == 0)
		cerr << "warning: `" << inFile << "': "
			"contains no usable sequence\n";

	if (opt::rank <= 0 && count == 0 && seqCollection->empty()) {
		assert(!opt::colourSpace);
		seqCollection->setColourSpace(false);
	}
}
/** The number of k-mer that have been eroded. */
size_t g_numEroded;

#if _SQL
std::vector<size_t> tempCounter(16,0);
InsOrderedMap<std::string,int> tempStatMap;

void addToDb(const std::string& key, const int& value)
{
	tempStatMap.push_back(key, value);
}
#endif

};
