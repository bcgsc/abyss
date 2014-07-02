#include "AssemblyAlgorithms.h"
#include "Assembly/Options.h"
#include "Common/Options.h"
#include "FastaReader.h"
#include "FastaWriter.h"
#include "Histogram.h"
#include "IOUtil.h"
#include "Log.h"
#include "SequenceCollection.h"
#include "StringUtil.h"
#include "Timer.h"
#include <algorithm>
#include <cctype>
#include <climits> // for UINT_MAX
#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;

namespace AssemblyAlgorithms
{
#if _SQL
vector<size_t> tempCounter(16,0);
InsOrderedMap<string,int> tempStatMap;

void addToDb(const string& key, const int& value)
{
	tempStatMap.push_back(key, value);
}
#endif

/** Return the kmer which are adjacent to this kmer. */
void generateSequencesFromExtension(const Kmer& currSeq,
		extDirection dir, SeqExt extension, vector<Kmer>& outseqs)
{
	vector<Kmer> extensions;
	Kmer extSeq(currSeq);
	extSeq.shift(dir);

	// Check for the existance of the 4 possible extensions
	for (unsigned i = 0; i < NUM_BASES; i++) {
		// Does this sequence have an extension?
		if(extension.checkBase(i))
		{
			extSeq.setLastBase(dir, i);
			outseqs.push_back(extSeq);
		}
	}
}

/** Load k-mer with coverage data.
 * @return the number of k-mer loaded
 */
static size_t loadKmer(ISequenceCollection& g, FastaReader& in)
{
	assert(opt::rank == -1);
	size_t count = 0;
	for (FastaRecord rec; in >> rec;) {
		assert(rec.seq.size() == opt::kmerSize);
		istringstream iss(rec.id);
		float coverage = 1;
		iss >> coverage;
		assert(iss);
		assert(iss.eof());
		g.add(Kmer(rec.seq), max(1, (int)ceilf(coverage)));

		if (++count % 1000000 == 0) {
			logger(1) << "Read " << count << " k-mer. ";
			g.printLoad();
		}
		g.pumpNetwork();
	}
	assert(in.eof());
	return count;
}

/** Load sequence data into the collection. */
void loadSequences(ISequenceCollection* seqCollection, string inFile)
{
	Timer timer("LoadSequences " + inFile);

	logger(0) << "Reading `" << inFile << "'...\n";

	if (inFile.find(".kmer") != string::npos) {
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
		if (opt::kmerSize > len) {
			count_small++;
			continue;
		}

		if (opt::rank <= 0
				&& count == 0 && seqCollection->empty()) {
			// Detect colour-space reads.
			bool colourSpace
				= seq.find_first_of("0123") != string::npos;
			seqCollection->setColourSpace(colourSpace);
			if (colourSpace)
				cout << "Colour-space assembly\n";
		}

		if (isalnum(seq[0])) {
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));
		}

		bool good = seq.find_first_not_of("ACGT0123") == string::npos;
		bool discarded = true;

		if (opt::ss && rec.id.size() > 2
				&& rec.id.substr(rec.id.size()-2) == "/1") {
			seq = reverseComplement(seq);
			count_reversed++;
		}

		for (unsigned i = 0; i < len - opt::kmerSize + 1; i++) {
			Sequence kmer(seq, i, opt::kmerSize);
			if (good || kmer.find_first_not_of("acgtACGT0123")
					== string::npos) {
				if (good || kmer.find_first_of("acgt") == string::npos)
					seqCollection->add(Kmer(kmer));
				else {
					transform(kmer.begin(), kmer.end(), kmer.begin(),
							::toupper);
					seqCollection->add(Kmer(kmer), 0);
				}
				discarded = false;
			}
		}
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
		cerr << "`" << inFile << "': "
			"reversed " << count_reversed << " reads\n";
	if (count_small > 0)
		cerr << "`" << inFile << "': "
			"discarded " << count_small << " reads "
			"shorter than " << opt::kmerSize << " bases\n";
#if _SQL
			tempCounter[10] += count_small;
#endif
	if (reader.unchaste() > 0)
		cerr << "`" << inFile << "': "
			"discarded " << reader.unchaste() << " unchaste reads\n";
	if (count_nonACGT > 0)
		cerr << "`" << inFile << "': "
			"discarded " << count_nonACGT << " reads "
			"containing non-ACGT characters\n";
#if _SQL
			tempCounter[11] += count_nonACGT;
#endif
	if (count_good == 0)
		cerr << "warning: `" << inFile << "': "
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
}

/** Generate the adjacency information for each sequence in the
 * collection. */
size_t generateAdjacency(ISequenceCollection* seqCollection)
{
	Timer timer("GenerateAdjacency");

	size_t count = 0;
	size_t numBasesSet = 0;
	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		if (++count % 1000000 == 0)
			logger(1) << "Finding adjacent k-mer: " << count << '\n';

		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			Kmer testSeq(iter->first);
			uint8_t adjBase = testSeq.shift(dir);
			for (unsigned i = 0; i < NUM_BASES; i++) {
				testSeq.setLastBase(dir, i);
				if (seqCollection->setBaseExtension(
							testSeq, !dir, adjBase))
					numBasesSet++;
			}
		}
		seqCollection->pumpNetwork();
	}

	if (numBasesSet > 0) {
		logger(0) << "Added " << numBasesSet << " edges.\n";
#if _SQL
		addToDb("EdgesGenerated", numBasesSet);
#endif
	}
	return numBasesSet;
}

/** Mark the specified vertex and its neighbours.
 * @return the number of marked edges
 */
static size_t markNeighbours(ISequenceCollection* g,
		const ISequenceCollection::value_type& u, extDirection sense)
{
	vector<Kmer> adj;
	generateSequencesFromExtension(u.first, sense,
			u.second.getExtension(sense), adj);
	for (vector<Kmer>::iterator v = adj.begin(); v != adj.end(); ++v)
		g->mark(*v, !sense);
	return adj.size();
}

/** Mark ambiguous branches and branches from palindromes for removal.
 * @return the number of branches marked
 */
size_t markAmbiguous(ISequenceCollection* g)
{
	Timer timer(__func__);
	size_t progress = 0;
	size_t countv = 0, counte = 0;
	for (ISequenceCollection::iterator it = g->begin();
			it != g->end(); ++it) {
		if (it->second.deleted())
			continue;

		if (++progress % 1000000 == 0)
			logger(1) << "Splitting: " << progress << '\n';

		if (!opt::ss && it->first.isPalindrome()) {
			countv += 2;
			g->mark(it->first);
			counte += markNeighbours(g, *it, SENSE);
		} else {
			for (extDirection sense = SENSE;
					sense <= ANTISENSE; ++sense) {
				if (it->second.getExtension(sense).isAmbiguous()
						|| (!opt::ss && it->first.isPalindrome(sense))) {
					countv++;
					g->mark(it->first, sense);
					counte += markNeighbours(g, *it, sense);
				}
			}
		}

		g->pumpNetwork();
	}
#if _SQL
	tempCounter[5] = countv;
	tempCounter[6] = counte;
#endif
	logger(0) << "Marked " << counte << " edges of " << countv
		<< " ambiguous vertices." << endl;

	return countv;
}

/** Remove the edges of marked and deleted vertices.
 * @return the number of branches removed
 */
size_t splitAmbiguous(ISequenceCollection* pSC)
{
	Timer timer(__func__);
	size_t count = 0;
	for (ISequenceCollection::iterator it = pSC->begin();
			it != pSC->end(); ++it) {
		if (!it->second.deleted())
			continue;
		for (extDirection sense = SENSE;
				sense <= ANTISENSE; ++sense) {
			if (it->second.marked(sense)) {
				removeExtensionsToSequence(pSC, *it, sense);
				count++;
			}
		}
		pSC->pumpNetwork();
	}
#if _SQL
	tempCounter[7] += count;
#endif
	logger(0) << "Split " << count << " ambigiuous branches.\n";
	return count;
}

/** Open the bubble file. */
void openBubbleFile(ofstream& out)
{
	if (opt::snpPath.empty())
		return;
	string path;
	if (opt::rank < 0) {
		path = opt::snpPath;
	} else {
		ostringstream s;
		s << "snp-" << opt::rank << ".fa";
		path = s.str();
	}
	out.open(path.c_str());
	assert_good(out, path);
}

/** Pop bubbles. */
size_t popBubbles(SequenceCollectionHash* seqCollection, ostream& out)
{
	Timer timer("PopBubbles");
	size_t numPopped = 0;

	// Set the cutoffs
	const unsigned maxNumBranches = 3;
	const unsigned maxLength = opt::bubbleLen - opt::kmerSize + 1;

	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		ExtensionRecord extRec = iter->second.extension();
		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			if (extRec.dir[dir].isAmbiguous()) {
				// Found a potential bubble, examine each branch
				bool stop = false;

				// Create the branch group
				BranchGroup branchGroup(dir, maxNumBranches,
						iter->first);
				initiateBranchGroup(branchGroup, iter->first,
						extRec.dir[dir]);

				// Iterate over the branches
				while(!stop)
				{
					size_t numBranches = branchGroup.size();
					for (unsigned j = 0; j < numBranches; ++j) {
						// Get the extensions of this branch
						ExtensionRecord extRec;
						int multiplicity = -1;

						const Kmer& lastKmer
							= branchGroup[j].back().first;
						bool success = seqCollection->getSeqData(
								lastKmer, extRec, multiplicity);
						assert(success);
						(void)success;
						processBranchGroupExtension(branchGroup, j,
								lastKmer, extRec, multiplicity,
								maxLength);
					}

					// At this point all branches should have the same
					// length or one will be a noext.
					branchGroup.updateStatus(maxLength);
					BranchGroupStatus status
						= branchGroup.getStatus();
					if (status == BGS_TOOLONG
							|| status == BGS_TOOMANYBRANCHES
							|| status == BGS_NOEXT) {
						stop = true;
					}
					else if(status == BGS_JOINED)
					{
						static unsigned snpID;
						writeBubble(out, branchGroup, ++snpID);
						assert(branchGroup.isAmbiguous(
									*seqCollection));
						collapseJoinedBranches(seqCollection,
								branchGroup);
						assert(!branchGroup.isAmbiguous(
									*seqCollection));
						numPopped++;
						stop = true;
					} else
						assert(status == BGS_ACTIVE);
				}
			}
		}
		seqCollection->pumpNetwork();
	}

	if (numPopped > 0)
		cout << "Removed " << numPopped << " bubbles.\n";
#if _SQL
	addToDb("totalErodedTips", tempCounter[0]);
	addToDb("totalPrunedTips", tempCounter[1]);
	addToDb("prunedRounds", tempCounter[2]);
	addToDb("totalLowCovCntg", tempCounter[3]);
	addToDb("totalLowCovKmer", tempCounter[4]);
	addToDb("totalSplitAmbg", tempCounter[7]);
	addToDb("poppedBubbles", numPopped);
	tempCounter.assign(8,0);
#endif
	return numPopped;
}

// Populate a branch group with the inital branches from a sequence
void initiateBranchGroup(BranchGroup& group, const Kmer& seq,
		const SeqExt& extension)
{
	vector<Kmer> extSeqs;
	generateSequencesFromExtension(seq, group.getDirection(),
			extension, extSeqs);
	assert(extSeqs.size() > 1);
	for (vector<Kmer>::iterator seqIter = extSeqs.begin();
			seqIter != extSeqs.end(); ++seqIter)
		group.addBranch(BranchRecord(group.getDirection()), *seqIter);
}

/** Process an a branch group extension. */
bool processBranchGroupExtension(BranchGroup& group,
		size_t branchIndex, const Kmer& seq,
		ExtensionRecord ext, int multiplicity,
		unsigned maxLength)
{
	BranchRecord& branch = group[branchIndex];
	branch.setData(make_pair(seq, KmerData(multiplicity, ext)));

	extDirection dir = group.getDirection();
	if (ext.dir[!dir].isAmbiguous()) {
		// Check that this fork is due to branches of our bubble
		// merging back together. If not, stop this bubble.
		if (branch.size() < 2) {
			group.setNoExtension();
			return false;
		}

		vector<Kmer> extKmer;
		generateSequencesFromExtension(seq, !dir,
				ext.dir[!dir], extKmer);
		assert(extKmer.size() > 1);
		for (vector<Kmer>::iterator it = extKmer.begin();
				it != extKmer.end(); ++it) {
			assert(branch.size() > 1);
			if (!group.exists(branch.size() - 2, *it)) {
				group.setNoExtension();
				return false;
			}
		}
		// Ignore the ambiguity.
		ext.dir[!dir].clear();
	}

	if (ext.dir[dir].isAmbiguous()) {
		// Create a new branch to follow the fork.
		vector<Kmer> extKmer;
		generateSequencesFromExtension(seq, dir,
				ext.dir[dir], extKmer);
		assert(extKmer.size() > 1);
		BranchRecord original = branch;
		vector<Kmer>::iterator it = extKmer.begin();
		branch.push_back(make_pair(*it++, KmerData()));
		for (; it != extKmer.end(); ++it)
			group.addBranch(original, *it);
		return group.isExtendable();
	}

	Kmer nextKmer = seq;
	if (processLinearExtensionForBranch(branch,
			nextKmer, ext, multiplicity,
			maxLength, false))
		branch.push_back(make_pair(nextKmer, KmerData()));
	else
		group.setNoExtension();
	return group.isExtendable();
}

/** Write a bubble to the specified file. */
void writeBubble(ostream& out, const BranchGroup& group, unsigned id)
{
	if (opt::snpPath.empty())
		return;

	char allele = 'A';
	for (BranchGroup::const_iterator it = group.begin();
			it != group.end(); ++it) {
		const BranchRecord& currBranch = *it;
		Sequence contig(currBranch);
		out << '>' << id << allele++ << ' '
			<< contig.length() << ' '
			<< currBranch.calculateBranchMultiplicity() << '\n'
			<< contig.c_str() << '\n';
	}
	assert(out.good());
}

/** Collapse a bubble to a single path. */
void collapseJoinedBranches(ISequenceCollection* collection,
		BranchGroup& group)
{
	const BranchRecord& best = group[0];
	logger(5) << "Popping " << best.size() << ' '
		<< best.front().first << '\n';

	// Add the k-mer from the dead branches.
	map<Kmer, KmerData> doomed;
	for (BranchGroup::const_iterator branchIt = group.begin() + 1;
			branchIt != group.end(); ++branchIt) {
		const BranchRecord& branch = *branchIt;
		for (BranchRecord::const_iterator it = branch.begin();
				it != branch.end(); ++it)
			doomed.insert(*it);
	}

	// Remove the k-mer that are in the good branch.
	for (BranchRecord::const_iterator it = best.begin();
			it != best.end(); ++it)
		doomed.erase(it->first);

	// Remove the dead k-mer from the assembly.
	for (map<Kmer, KmerData>::const_iterator it = doomed.begin();
			it != doomed.end(); ++it)
		removeSequenceAndExtensions(collection, *it);
}

/**
 * Remove a k-mer and update the extension records of the k-mer that
 * extend to it.
 */
void removeSequenceAndExtensions(ISequenceCollection* seqCollection,
		const ISequenceCollection::value_type& seq)
{
	// This removes the reverse complement as well
	seqCollection->remove(seq.first);
	removeExtensionsToSequence(seqCollection, seq, SENSE);
	removeExtensionsToSequence(seqCollection, seq, ANTISENSE);
}

/** Remove all the extensions to this sequence. */
void removeExtensionsToSequence(ISequenceCollection* seqCollection,
		const ISequenceCollection::value_type& seq, extDirection dir)
{
	SeqExt extension(seq.second.getExtension(dir));
	Kmer testSeq(seq.first);
	uint8_t extBase = testSeq.shift(dir);
	for (unsigned i = 0; i < NUM_BASES; i++) {
		if (extension.checkBase(i)) {
			testSeq.setLastBase(dir, i);
			seqCollection->removeExtension(testSeq, !dir, extBase);
		}
	}
}

/** The number of k-mer that have been eroded. */
static size_t g_numEroded;

/** Return the number of k-mer that have been eroded. */
size_t getNumEroded()
{
	size_t numEroded = g_numEroded;
	g_numEroded = 0;
#if _SQL
	tempCounter[0] += numEroded;
#endif
	logger(0) << "Eroded " << numEroded << " tips.\n";
	return numEroded;
}

/** Consider the specified k-mer for erosion.
 * @return the number of k-mer eroded, zero or one
 */
size_t erode(ISequenceCollection* c,
		const ISequenceCollection::value_type& seq)
{
	if (seq.second.deleted())
		return 0;
	extDirection dir;
	SeqContiguity contiguity = checkSeqContiguity(seq, dir);
	if (contiguity == SC_CONTIGUOUS)
		return 0;

	const KmerData& data = seq.second;
	if (data.getMultiplicity() < opt::erode
			|| data.getMultiplicity(SENSE) < opt::erodeStrand
			|| data.getMultiplicity(ANTISENSE) < opt::erodeStrand) {
		removeSequenceAndExtensions(c, seq);
		g_numEroded++;
		return 1;
	} else
		return 0;
}

/** The given sequence has changed. */
static void erosionObserver(ISequenceCollection* c,
		const ISequenceCollection::value_type& seq)
{
	erode(c, seq);
}

//
// Erode data off the ends of the graph, one by one
//
size_t erodeEnds(ISequenceCollection* seqCollection)
{
	Timer erodeEndsTimer("Erode");
	assert(g_numEroded == 0);
	seqCollection->attach(erosionObserver);

	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		erode(seqCollection, *iter);
		seqCollection->pumpNetwork();
	}

	seqCollection->detach(erosionObserver);
	return getNumEroded();
}

static size_t trimSequences(SequenceCollectionHash* seqCollection,
		unsigned maxBranchCull);

/** Trimming driver function */
void performTrim(SequenceCollectionHash* seqCollection)
{
	if (opt::trimLen == 0)
		return;
	unsigned rounds = 0;
	size_t total = 0;
	for (unsigned trim = 1; trim < opt::trimLen; trim *= 2) {
		rounds++;
		total += trimSequences(seqCollection, trim);
	}
	size_t count;
	while ((count = trimSequences(seqCollection, opt::trimLen)) > 0) {
		rounds++;
		total += count;
	}
	cout << "Pruned " << total << " tips in "
		<< rounds << " rounds.\n";
#if _SQL
	tempCounter[1] += total;
	tempCounter[2] = rounds;
#endif
}

/** Return the adjacency of this sequence.
 * @param considerMarks when true, treat a marked vertex as having
 * no edges
 */
SeqContiguity checkSeqContiguity(
		const ISequenceCollection::value_type& seq,
		extDirection& outDir, bool considerMarks)
{
	assert(!seq.second.deleted());
	bool child = seq.second.hasExtension(SENSE)
		&& !(considerMarks && seq.second.marked(SENSE));
	bool parent = seq.second.hasExtension(ANTISENSE)
		&& !(considerMarks && seq.second.marked(ANTISENSE));
	if(!child && !parent)
	{
		//this sequence is completely isolated
		return SC_ISLAND;
	}
	else if(!child)
	{
		outDir = ANTISENSE;
		return SC_ENDPOINT;
	}
	else if(!parent)
	{
		outDir = SENSE;
		return SC_ENDPOINT;
	}
	else
	{
		// sequence is contiguous
		return SC_CONTIGUOUS;
	}
}

/** Prune tips shorter than maxBranchCull. */
static size_t trimSequences(SequenceCollectionHash* seqCollection,
		unsigned maxBranchCull)
{
	Timer timer("TrimSequences");
	cout << "Pruning tips shorter than "
		<< maxBranchCull << " bp...\n";
	size_t numBranchesRemoved = 0;

	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		extDirection dir;
		// dir will be set to the trimming direction if the sequence
		// can be trimmed.
		SeqContiguity status = checkSeqContiguity(*iter, dir);

		if (status == SC_CONTIGUOUS)
			continue;
		else if(status == SC_ISLAND)
		{
			// remove this sequence, it has no extensions
			seqCollection->mark(iter->first);
			numBranchesRemoved++;
			continue;
		}

		BranchRecord currBranch(dir);
		Kmer currSeq = iter->first;
		while(currBranch.isActive())
		{
			ExtensionRecord extRec;
			int multiplicity = -1;
			bool success = seqCollection->getSeqData(
					currSeq, extRec, multiplicity);
			assert(success);
			(void)success;
			processLinearExtensionForBranch(currBranch,
					currSeq, extRec, multiplicity, maxBranchCull);
		}

		// The branch has ended check it for removal, returns true if
		// it was removed.
		if(processTerminatedBranchTrim(seqCollection, currBranch))
		{
			numBranchesRemoved++;
		}
		seqCollection->pumpNetwork();
	}

	size_t numSweeped = removeMarked(seqCollection);

	if (numBranchesRemoved > 0)
		logger(0) << "Pruned " << numSweeped << " k-mer in "
			<< numBranchesRemoved << " tips.\n";
	return numBranchesRemoved;
}

/** Extend this branch. */
bool extendBranch(BranchRecord& branch, Kmer& kmer, SeqExt ext)
{
	if (!ext.hasExtension()) {
		branch.terminate(BS_NOEXT);
		return false;
	} else if (ext.isAmbiguous()) {
		branch.terminate(BS_AMBI_SAME);
		return false;
	} else {
		vector<Kmer> adj;
		generateSequencesFromExtension(kmer, branch.getDirection(),
				ext, adj);
		assert(adj.size() == 1);
		kmer = adj.front();
		return true;
	}
}

/**
 * Process the extension for this branch for the trimming algorithm
 * CurrSeq is the current sequence being inspected (the next member to
 * be added to the branch). The extension record is the extensions of
 * that sequence and multiplicity is the number of times that kmer
 * appears in the data set. After processing currSeq is unchanged if
 * the branch is no longer active or else it is the generated
 * extension. If the parameter addKmer is true, add the k-mer to the
 * branch.
 */
bool processLinearExtensionForBranch(BranchRecord& branch,
		Kmer& currSeq, ExtensionRecord extensions, int multiplicity,
		unsigned maxLength, bool addKmer)
{
	/** Stop contig assembly at palindromes. */
	const bool stopAtPalindromes = !opt::ss && maxLength == UINT_MAX;

	extDirection dir = branch.getDirection();
	if (branch.isTooLong(maxLength)) {
		// Too long.
		branch.terminate(BS_TOO_LONG);
		return false;
	} else if (extensions.dir[!dir].isAmbiguous()) {
		// Ambiguous.
		branch.terminate(BS_AMBI_OPP);
		return false;
	} else if (stopAtPalindromes && currSeq.isPalindrome()) {
		// Palindrome.
		branch.terminate(BS_AMBI_SAME);
		return false;
	}

	if (addKmer)
		branch.push_back(make_pair(currSeq,
					KmerData(multiplicity, extensions)));

	if (branch.isTooLong(maxLength)) {
		// Too long.
		branch.terminate(BS_TOO_LONG);
		return false;
	} else if (stopAtPalindromes && currSeq.isPalindrome(dir)) {
		// Palindrome.
		branch.terminate(BS_AMBI_SAME);
		return false;
	}

	return extendBranch(branch, currSeq, extensions.dir[dir]);
}

/** Trim the specified branch if it meets trimming criteria.
 * @return true if the specified branch was trimmed
 */
bool processTerminatedBranchTrim(ISequenceCollection* seqCollection,
		BranchRecord& branch)
{
	assert(!branch.isActive());
	assert(!branch.empty());
	if (branch.getState() == BS_NOEXT
			|| branch.getState() == BS_AMBI_OPP) {
		logger(5) << "Pruning " << branch.size() << ' '
			<< branch.front().first << '\n';
		for (BranchRecord::iterator it = branch.begin();
				it != branch.end(); ++it)
			seqCollection->mark(it->first);
		return true;
	} else
		return false;
}

/** Remove all marked k-mer.
 * @return the number of removed k-mer
 */
size_t removeMarked(ISequenceCollection* pSC)
{
	Timer timer(__func__);
	size_t count = 0;
	for (ISequenceCollection::iterator it = pSC->begin();
			it != pSC->end(); ++it) {
		if (it->second.deleted())
			continue;
		if (it->second.marked()) {
			removeSequenceAndExtensions(pSC, *it);
			count++;
		}
		pSC->pumpNetwork();
	}
	if (count > 0)
		logger(1) << "Removed " << count << " marked k-mer.\n";
	return count;
}

/** Assemble a contig.
 * @return the number of k-mer below the coverage threshold
 */
size_t assembleContig(
		ISequenceCollection* seqCollection, FastaWriter* writer,
		BranchRecord& branch, unsigned id)
{
	assert(!branch.isActive());
	assert(branch.getState() == BS_NOEXT
			|| branch.getState() == BS_AMBI_SAME
			|| branch.getState() == BS_AMBI_OPP);

	// Assemble the contig.
	Sequence contig(branch);

	size_t kmerCount = branch.calculateBranchMultiplicity();
	if (writer != NULL)
		writer->WriteSequence(contig, id, kmerCount);

	// Remove low-coverage contigs.
	float coverage = (float)kmerCount / branch.size();
	if (opt::coverage > 0 && coverage < opt::coverage) {
		for (BranchRecord::iterator it = branch.begin();
				it != branch.end(); ++it)
			seqCollection->remove(it->first);
		return branch.size();
	}
	return 0;
}

/** Assemble contigs.
 * @return the number of contigs assembled
 */
size_t assemble(SequenceCollectionHash* seqCollection,
		FastaWriter* fileWriter)
{
	Timer timer("Assemble");

	size_t kmerCount = 0;
	unsigned contigID = 0;
	size_t assembledKmer = 0;
	size_t lowCoverageKmer = 0;
	size_t lowCoverageContigs = 0;

	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;
		kmerCount++;

		extDirection dir;
		SeqContiguity status = checkSeqContiguity(*iter, dir, true);
		if (status == SC_CONTIGUOUS)
			continue;
		else if(status == SC_ISLAND)
		{
			BranchRecord currBranch(SENSE);
			currBranch.push_back(*iter);
			currBranch.terminate(BS_NOEXT);
			size_t removed = assembleContig(seqCollection,
					fileWriter, currBranch, contigID++);
			assembledKmer += currBranch.size();
			if (removed > 0) {
				lowCoverageContigs++;
				lowCoverageKmer += removed;
			}
			continue;
		}
		assert(status == SC_ENDPOINT);

		BranchRecord currBranch(dir);
		currBranch.push_back(*iter);
		Kmer currSeq = iter->first;
		extendBranch(currBranch, currSeq,
				iter->second.getExtension(dir));
		assert(currBranch.isActive());
		while(currBranch.isActive())
		{
			ExtensionRecord extRec;
			int multiplicity = -1;
			bool success = seqCollection->getSeqData(
					currSeq, extRec, multiplicity);
			assert(success);
			(void)success;
			processLinearExtensionForBranch(currBranch,
					currSeq, extRec, multiplicity, UINT_MAX);
		}

		if ((opt::ss && currBranch.getDirection() == SENSE)
				|| (!opt::ss && currBranch.isCanonical())) {
			size_t removed = assembleContig(seqCollection,
					fileWriter, currBranch, contigID++);
			assembledKmer += currBranch.size();
			if (removed > 0) {
				lowCoverageContigs++;
				lowCoverageKmer += removed;
			}
		}

		seqCollection->pumpNetwork();
	}

	if (opt::coverage > 0) {
		cout << "Found " << assembledKmer << " k-mer in " << contigID
			<< " contigs before removing low-coverage contigs.\n"
			"Removed " << lowCoverageKmer << " k-mer in "
				<< lowCoverageContigs << " low-coverage contigs.\n";
#if _SQL
		tempCounter[3] += lowCoverageContigs;
		tempCounter[4] += lowCoverageKmer;
#endif
	} else {
		assert(assembledKmer <= kmerCount);
		size_t circularKmer = kmerCount - assembledKmer;
		if (circularKmer > 0)
			cout << "Left " << circularKmer
				<< " unassembled k-mer in circular contigs.\n";
		cout << "Assembled " << assembledKmer << " k-mer in "
			<< contigID << " contigs.\n";
#if _SQL
		addToDb("finalAmbgVertices", tempCounter[5]);
		addToDb("finalAmbgEdges", tempCounter[6]);
		tempCounter.assign(8,0);
		addToDb("unassembledCircularCntg", circularKmer);
		addToDb("assembledKmerNum", assembledKmer);
		addToDb("assembledCntg", contigID);
#endif
	}
	return contigID;
}

/** Return the k-mer coverage histogram. */
Histogram coverageHistogram(const ISequenceCollection& c)
{
	Histogram h;
	for (ISequenceCollection::const_iterator it = c.begin();
			it != c.end(); ++it) {
		if (it->second.deleted())
			continue;
		h.insert(it->second.getMultiplicity());
	}
	return h;
}

/** Calculate a k-mer coverage threshold from the given k-mer coverage
 * histogram. */
static float calculateCoverageThreshold(const Histogram& h)
{
	float cov = h.firstLocalMinimum();
	if (opt::rank <= 0) {
		if (cov == 0)
			cout << "Unable to determine minimum k-mer coverage\n";
		else
			cout << "Minimum k-mer coverage is " << cov << endl;
	}

	for (unsigned iteration = 0; iteration < 100; iteration++) {
		Histogram trimmed = h.trimLow((unsigned)roundf(cov));
		if (opt::rank <= 0)
			logger(1) << "Coverage: " << cov << "\t"
				"Reconstruction: " << trimmed.size() << endl;

		unsigned median = trimmed.median();
		float cov1 = sqrt(median);
		if (cov1 == cov) {
			// The coverage threshold has converged.
			if (opt::rank <= 0)
				cout << "Using a coverage threshold of "
					<< (unsigned)roundf(cov) << "...\n"
					"The median k-mer coverage is " << median << "\n"
					"The reconstruction is " << trimmed.size()
					<< endl;
#if _SQL
			addToDb("coverageThreshold", (unsigned)roundf(cov));
			addToDb("medianKcoverage", median);
			addToDb("restruction", trimmed.size());
#endif
			return cov;
		}
		cov = cov1;
	}
	if (opt::rank <= 0)
		cerr << "warning: coverage threshold did not converge"
			<< endl;
	return 0;
}

/** Set the coverage-related parameters e and c from the given k-mer
 * coverage histogram. */
void setCoverageParameters(const Histogram& h)
{
	if (!opt::coverageHistPath.empty() && opt::rank <= 0) {
		ofstream histFile(opt::coverageHistPath.c_str());
		assert_good(histFile, opt::coverageHistPath);
		histFile << h;
		assert(histFile.good());
	}

	float minCov = calculateCoverageThreshold(h);
	if (opt::rank <= 0) {
		if (minCov == 0)
			cout << "Unable to determine the "
				"k-mer coverage threshold" << endl;
		else
			cout << "The k-mer coverage threshold is " << minCov
				<< endl;
	}
	if (minCov < 2)
		minCov = 2;

	if ((int)opt::erode < 0) {
		opt::erode = (unsigned)roundf(minCov);
		if (opt::rank <= 0)
			cout << "Setting parameter e (erode) to "
				<< opt::erode << endl;
	}
	if ((int)opt::erodeStrand < 0) {
		opt::erodeStrand = minCov <= 2 ? 0 : 1;
		if (opt::rank <= 0)
			cout << "Setting parameter E (erodeStrand) to "
				<< opt::erodeStrand << endl;
	}
	if (opt::coverage < 0) {
		opt::coverage = minCov;
		if (opt::rank <= 0)
			cout << "Setting parameter c (coverage) to "
				<< opt::coverage << endl;
	}
}

};
