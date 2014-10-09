#ifndef PAIREDDBG_LOADALGORITHM_H
#define PAIREDDBG_LOADALGORITHM_H 1

#include "DataLayer/FastaReader.h"

namespace AssemblyAlgorithms {

/** Load k-mer with coverage data.
 * @return the number of k-mer loaded
 */
static inline
size_t loadKmer(ISequenceCollection& g, FastaReader& in)
{
	typedef SequenceCollectionHash Graph;
	typedef graph_traits<Graph>::vertex_descriptor V;

	assert(opt::rank == -1);
	size_t count = 0;
	for (FastaRecord rec; in >> rec;) {
		assert(rec.seq.size() == V::length());
		std::istringstream iss(rec.id);
		float coverage = 1;
		iss >> coverage;
		assert(iss);
		assert(iss.eof());

		g.add(V(rec.seq), std::max(1, (int)ceilf(coverage)));

		if (++count % 1000000 == 0) {
			logger(1) << "Read " << count << " k-mer. ";
			g.printLoad();
		}
		g.pumpNetwork();
	}
	assert(in.eof());
	return count;
}

static inline
bool loadSequence(ISequenceCollection* seqCollection, Sequence& seq)
{
	size_t len = seq.length();

	if (isalnum(seq[0])) {
		if (opt::colourSpace)
			assert(isdigit(seq[0]));
		else
			assert(isalpha(seq[0]));
	}

	bool good = seq.find_first_not_of("ACGT0123") == std::string::npos;
	bool discarded = true;

	for (unsigned i = 0; i < len - KmerPair::length() + 1; ++i) {
		Sequence kmer;
		int gapSize = KmerPair::length() - 2 * Kmer::length();
		if (gapSize > 0)
		{
			std::stringstream sseed;
			sseed << Sequence(seq, i, Kmer::length()) << Sequence(seq, i + Kmer::length() + gapSize, Kmer::length());
			kmer = sseed.str();
			if (sseed.str().size() != 2 * Kmer::length()) {
				std::cout << Sequence(seq, i, Kmer::length() + gapSize)
					<< '\n'
					<< Sequence(seq, i, Kmer::length()) << '\n'
					<< Sequence(seq, i + Kmer::length() + gapSize, Kmer::length())
					<< '\n';
				exit(EXIT_FAILURE);
			}
		}
		// if delta=0
		else
			kmer = Sequence(seq, i, KmerPair::length());

		if (good || kmer.find_first_not_of("acgtACGT0123") == std::string::npos)
		{
			if (good || kmer.find_first_of("acgt") == std::string::npos)
			{
				Sequence buf1(kmer, 0, Kmer::length());
				Sequence buf2(kmer, Kmer::length(), Kmer::length());
				seqCollection->add(KmerPair(buf1,buf2));
			}
			else
			{
				transform(kmer.begin(), kmer.end(), kmer.begin(),::toupper);
				Sequence buf1(kmer, 0, Kmer::length());
				Sequence buf2(kmer, Kmer::length(), Kmer::length());
				seqCollection->add(KmerPair(buf1,buf2), 0);
			}
			discarded = false;
		}
	}

	return discarded;

}

/** Load sequence data into the collection. */
static inline
void loadSequences(ISequenceCollection* seqCollection, std::string inFile)
{
	Timer timer("LoadSequences " + inFile);

	logger(0) << "Reading `" << inFile << "'...\n";

	if (inFile.find(".kmer") != std::string::npos) {
		if (opt::rank <= 0)
			seqCollection->setColourSpace(false);
		seqCollection->load(inFile.c_str());
		return;
	}

	size_t count = 0, count_good = 0,
			 count_small = 0, count_nonACGT = 0,
			 count_reversed = 0;
	int fastaFlags = opt::maskCov ?  FastaReader::NO_FOLD_CASE :
			FastaReader::FOLD_CASE;
	FastaReader reader(inFile.c_str(), fastaFlags);
	if (endsWith(inFile, ".jf") || endsWith(inFile, ".jfq")) {
		// Load k-mer with coverage data.
		count = loadKmer(*seqCollection, reader);
		count_good = count;
	} else
	for (FastaRecord rec; reader >> rec;) {

		Sequence seq = rec.seq;
		size_t len = seq.length();

		if (KmerPair::length() > len) {
			count_small++;
			continue;
		}

		if (opt::rank <= 0
				&& count == 0 && seqCollection->empty()) {
			// Detect colour-space reads.
			bool colourSpace
				= seq.find_first_of("0123") != std::string::npos;
			seqCollection->setColourSpace(colourSpace);
			if (colourSpace)
				std::cout << "Colour-space assembly\n";
		}

		if (opt::ss && rec.id.size() > 2
				&& rec.id.substr(rec.id.size()-2) == "/1") {
			seq = reverseComplement(seq);
			count_reversed++;
		}

		bool discarded = loadSequence(seqCollection, seq);

		if (discarded)
			count_nonACGT++;
		else
			count_good++;

		if (++count % 100000 == 0) {
			logger(1) << "Read " << count << " reads. ";
			seqCollection->printLoad();
		}

		seqCollection->pumpNetwork();
	}
	assert(reader.eof());

	logger(1) << "Read " << count << " reads. ";
	seqCollection->printLoad();

	if (count_reversed > 0)
		std::cerr << "`" << inFile << "': "
			"reversed " << count_reversed << " reads\n";
	if (count_small > 0)
		std::cerr << "`" << inFile << "': "
			"discarded " << count_small << " reads "
			"shorter than " << KmerPair::length() << " bases\n";
	if (reader.unchaste() > 0)
		std::cerr << "`" << inFile << "': "
			"discarded " << reader.unchaste() << " unchaste reads\n";
	if (count_nonACGT > 0)
		std::cerr << "`" << inFile << "': "
			"discarded " << count_nonACGT << " reads "
			"containing non-ACGT characters\n";
#if _SQL
			tempCounter[0] += count_reversed;
			tempCounter[1] += (count_small + reader.unchaste() + count_nonACGT);
#endif
	if (count_good == 0)
		std::cerr << "warning: `" << inFile << "': "
			"contains no usable sequence\n";

	if (opt::rank <= 0 && count == 0 && seqCollection->empty()) {
		/* The master process did not load any data, which means that
		 * it hasn't told the slave processes whether this assembly is
		 * in colour-space. Rather than fail right now, assume that
		 * the assembly is not colour space. If the assumption is
		 * incorrect, the assembly will fail pretty quickly as soon as
		 * one of the slave processes sees a colour-space read.
		 */
		assert(!opt::colourSpace);
		seqCollection->setColourSpace(false);
	}
#if _SQL
	addToDb("reversedReads", tempCounter[0]);
	addToDb("totalDiscardedReads", tempCounter[1]);
	tempCounter.assign(2,0);
#endif
}

} // namespace AssemblyAlgorithms

#endif
