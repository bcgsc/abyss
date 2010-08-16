#include "AssemblyAlgorithms.h"
#include "Assembly/Options.h"
#include "Common/Options.h"
#include "FastaReader.h"
#include "FastaWriter.h"
#include "Histogram.h"
#include "ISequenceCollection.h"
#include "Log.h"
#include "Timer.h"
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring> // for strerror
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

namespace AssemblyAlgorithms
{

/** Return the kmer which are adjacent to this kmer. */
void generateSequencesFromExtension(const Kmer& currSeq,
		extDirection dir, SeqExt extension, vector<Kmer>& outseqs)
{
	vector<Kmer> extensions;
	Kmer extSeq(currSeq);
	extSeq.shift(dir);

	// Check for the existance of the 4 possible extensions
	for (int i  = 0; i < NUM_BASES; i++) {
		// Does this sequence have an extension?
		if(extension.checkBase(i))
		{
			extSeq.setLastBase(dir, i);
			outseqs.push_back(extSeq);
		}
	}
}

/** Load sequence data into the collection. */
void loadSequences(ISequenceCollection* seqCollection, string inFile)
{
	Timer timer("LoadSequences " + inFile);

	PrintDebug(0, "Reading `%s'\n", inFile.c_str());

	if (inFile.find(".kmer") != string::npos) {
		if (opt::rank <= 0)
			seqCollection->setColourSpace(false);
		seqCollection->load(inFile.c_str());
		return;
	}

	unsigned count = 0, count_good = 0,
			 count_small = 0, count_nonACGT = 0;
	FastaReader reader(inFile.c_str(), FastaReader::FOLD_CASE);
	for (Sequence seq; reader >> seq;) {
		int len = seq.length();
		if (opt::kmerSize > len) {
			count_small++;
			continue;
		}

		if (opt::rank <= 0
				&& count == 0 && seqCollection->count() == 0) {
			// Detect colour-space reads.
			bool colourSpace
				= seq.find_first_of("0123") != string::npos;
			seqCollection->setColourSpace(colourSpace);
			if (colourSpace)
				puts("Colour-space assembly");
		}

		if (isalnum(seq[0])) {
			if (opt::colourSpace)
				assert(isdigit(seq[0]));
			else
				assert(isalpha(seq[0]));
		}

		bool good = seq.find_first_not_of("ACGT0123") == string::npos;
		bool discarded = true;
		for (int i = 0; i < len - opt::kmerSize + 1; i++) {
			Sequence kmer(seq, i, opt::kmerSize);
			if (good || kmer.find_first_not_of("ACGT0123")
					== string::npos) {
				seqCollection->add(Kmer(kmer));
				discarded = false;
			}
		}
		if (discarded)
			count_nonACGT++;
		else
			count_good++;

		if (++count % 100000 == 0) {
			PrintDebug(1, "Read %u reads. ", count);
			seqCollection->printLoad();
		}
		seqCollection->pumpNetwork();
	}
	assert(reader.eof());

	PrintDebug(1, "Read %u reads. ", count);
	seqCollection->printLoad();

	if (count_small > 0)
		fprintf(stderr, "warning: discarded %u reads "
				"shorter than %u bases\n",
				count_small, opt::kmerSize);
	if (reader.unchaste() > 0)
		cerr << "warning: discarded " << reader.unchaste()
			<< " unchaste reads" << endl;
	if (count_nonACGT > 0)
		fprintf(stderr, "warning: discarded %u reads "
				"containing non-ACGT characters\n", count_nonACGT);

	if (count_good == 0)
		fprintf(stderr, "warning: `%s' contains no usable sequence\n",
				inFile.c_str());
}

//
// Generate the adjacency information for each sequence in the collection
//
void generateAdjacency(ISequenceCollection* seqCollection)
{
	Timer timer("GenerateAdjacency");

	unsigned count = 0;
	unsigned numBasesSet = 0;
	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		if (++count % 1000000 == 0)
			PrintDebug(1, "Generating adjacency: %u k-mer\n", count);

		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			Kmer testSeq(iter->first);
			uint8_t adjBase = testSeq.shift(dir);
			for (int i = 0; i < NUM_BASES; i++) {
				testSeq.setLastBase(dir, i);
				if (seqCollection->setBaseExtension(
							testSeq, !dir, adjBase))
					numBasesSet++;
			}
		}
		seqCollection->pumpNetwork();
	}

	if (numBasesSet > 0)
		PrintDebug(0, "Generated %u edges\n", numBasesSet);
}

/** Mark the specified vertex and its neighbours.
 * @return the number of marked edges
 */
static unsigned markVertexAndNeighbours(ISequenceCollection* g,
		const PackedSeq& u, extDirection sense)
{
	g->mark(u.first, sense);
	vector<Kmer> adj;
	generateSequencesFromExtension(u.first, sense,
			u.second.getExtension(sense), adj);
	assert(!adj.empty());
	for (vector<Kmer>::iterator v = adj.begin();
			v != adj.end(); ++v)
		g->mark(*v, !sense);
	return adj.size();
}

/** Mark ambiguous branches and branches from palindromes for removal.
 * @return the number of branches marked
 */
unsigned markAmbiguous(ISequenceCollection* g)
{
	Timer timer(__func__);
	unsigned progress = 0;
	unsigned countv = 0, counte = 0;
	for (ISequenceCollection::iterator it = g->begin();
			it != g->end(); ++it) {
		if (it->second.deleted())
			continue;

		if (++progress % 1000000 == 0)
			PrintDebug(1, "Splitting: %u k-mer\n", progress);

		if (it->first.isPalindrome()) {
			countv += 2;
			counte += markVertexAndNeighbours(g, *it, SENSE);
			counte += markVertexAndNeighbours(g, *it, ANTISENSE);
		} else {
			for (extDirection sense = SENSE;
					sense <= ANTISENSE; ++sense) {
				if (it->second.getExtension(sense).isAmbiguous()
						|| it->first.isPalindrome(sense)) {
					countv++;
					counte += markVertexAndNeighbours(g, *it, sense);
				}
			}
		}

		g->pumpNetwork();
	}
	logger(0) << "Marked " << counte << " edges of " << countv
		<< " ambiguous vertices." << endl;
	return countv;
}

/** Remove the edges of marked and deleted vertices.
 * @return the number of branches removed
 */
unsigned splitAmbiguous(ISequenceCollection* pSC)
{
	Timer timer(__func__);
	unsigned count = 0;
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
	PrintDebug(0, "Split %u ambiguous branches\n", count);
	return count;
}

static void assert_open(/*const*/ ofstream& f, const string& p)
{
	if (f.is_open())
		return;
	cerr << p << ": " << strerror(errno) << endl;
	exit(EXIT_FAILURE);
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
	assert_open(out, path);
}

int popBubbles(ISequenceCollection* seqCollection, ostream& out)
{
	Timer timer("PopBubbles");
	int numPopped = 0;

	// Set the cutoffs
	const unsigned int maxNumBranches = 3;

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
						extRec.dir[dir],
						opt::bubbleLen - opt::kmerSize + 1);

				// Iterate over the branches
				while(!stop)
				{
					size_t numBranches = branchGroup.size();
					for(unsigned int j = 0; j < numBranches; ++j)
					{						
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
								lastKmer, extRec, multiplicity);
					}

					// At this point all branches should have the same length or one will be a noext
					branchGroup.updateStatus();
					BranchGroupStatus status = branchGroup.getStatus();
					if (status == BGS_TOOLONG
							|| status == BGS_TOOMANYBRANCHES
							|| status == BGS_NOEXT) {
						stop = true;
					}
					else if(status == BGS_JOINED)
					{
						static unsigned snpID;
						writeBubble(out, branchGroup, ++snpID);
						collapseJoinedBranches(seqCollection,
								branchGroup);
						numPopped++;
						stop = true;
					}
					else
					{										
						// the branch is still active, continue
						assert(status == BGS_ACTIVE);
					}
				}
			}
		}
		seqCollection->pumpNetwork();
	}

	if (numPopped > 0)
		printf("Removed %u bubbles\n", numPopped);
	return numPopped;
}

// Populate a branch group with the inital branches from a sequence
void initiateBranchGroup(BranchGroup& group, const Kmer& seq,
		const SeqExt& extension, size_t maxBubbleSize)
{
	vector<Kmer> extSeqs;
	generateSequencesFromExtension(seq, group.getDirection(),
			extension, extSeqs);
	assert(extSeqs.size() > 1);
	for (vector<Kmer>::iterator seqIter = extSeqs.begin();
			seqIter != extSeqs.end(); ++seqIter) {
		BranchRecord newBranch(group.getDirection(), maxBubbleSize);
		group.addBranch(newBranch, *seqIter);
	}
}

/** Process an a branch group extension. */
bool processBranchGroupExtension(BranchGroup& group,
		size_t branchIndex, const Kmer& seq,
		ExtensionRecord ext, int multiplicity)
{
	BranchRecord& branch = group[branchIndex];
	branch.setData(make_pair(seq, KmerData(multiplicity, ext)));

	extDirection dir = group.getDirection();
	if (ext.dir[!dir].isAmbiguous()) {
		// Check that this fork is due to branches of our bubble
		// merging back together. If not, stop this bubble.
		vector<Kmer> extKmer;
		generateSequencesFromExtension(seq, !dir,
				ext.dir[!dir], extKmer);
		assert(extKmer.size() > 1);
		for (vector<Kmer>::iterator it = extKmer.begin();
				it != extKmer.end(); ++it) {
			if (!group.exists(*it)) {
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
			nextKmer, ext, multiplicity, false))
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
			<< currBranch.getBranchMultiplicity() << '\n'
			<< contig.c_str() << '\n';
	}
	assert(out.good());
}

/** Collapse a bubble to a single path. */
void collapseJoinedBranches(ISequenceCollection* collection,
		BranchGroup& group)
{
	assert(group.isAmbiguous(collection));

	const BranchRecord& best = group[0];
	PrintDebug(5, "Popping %zu %s\n", best.size(),
				best.front().first.decode().c_str());

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
	assert(!group.isAmbiguous(collection));
}

/**
 * Remove a k-mer and update the extension records of the k-mer that
 * extend to it.
 */
void removeSequenceAndExtensions(ISequenceCollection* seqCollection,
		const PackedSeq& seq)
{
	// This removes the reverse complement as well
	seqCollection->remove(seq.first);
	removeExtensionsToSequence(seqCollection, seq, SENSE);
	removeExtensionsToSequence(seqCollection, seq, ANTISENSE);
}

/** Remove all the extensions to this sequence. */
void removeExtensionsToSequence(ISequenceCollection* seqCollection,
		const PackedSeq& seq, extDirection dir)
{
	SeqExt extension(seq.second.getExtension(dir));
	extDirection oppDir = oppositeDirection(dir);
	Kmer testSeq(seq.first);
	uint8_t extBase = testSeq.shift(dir);
	for (int i = 0; i < NUM_BASES; i++) {
		if (extension.checkBase(i)) {
			testSeq.setLastBase(dir, i);
			seqCollection->removeExtension(testSeq, oppDir, extBase);
		}
	}
}

/** The number of k-mer that have been eroded. */
static unsigned g_numEroded;

/** Return the number of k-mer that have been eroded. */
unsigned getNumEroded()
{
	unsigned numEroded = g_numEroded;
	g_numEroded = 0;
	PrintDebug(0, "Eroded %u tips\n", numEroded);
	return numEroded;
}

/** Consider the specified k-mer for erosion.
 * @return the number of k-mer eroded, zero or one
 */
unsigned erode(ISequenceCollection* c, const PackedSeq& seq)
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
		const PackedSeq& seq)
{
	erode(c, seq);
}

//
// Erode data off the ends of the graph, one by one
//
unsigned erodeEnds(ISequenceCollection* seqCollection)
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

/** Trimming driver function */
void performTrim(ISequenceCollection* seqCollection, int start)
{
	if (opt::trimLen == 0)
		return;
	unsigned rounds = 0, total = 0;
	for (int trim = start; trim < opt::trimLen; trim *= 2) {
		rounds++;
		total += trimSequences(seqCollection, trim);
	}
	unsigned count;
	while ((count = trimSequences(seqCollection, opt::trimLen)) > 0) {
		rounds++;
		total += count;
	}
	printf("Trimmed %u branches in %u rounds\n", total, rounds);
}

/** Return the adjacency of this sequence.
 * @param considerMarks when true, treat a marked vertex as having
 * no edges
 */
SeqContiguity checkSeqContiguity(const PackedSeq& seq,
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

//
// Trimming (error removal) function
//
int trimSequences(ISequenceCollection* seqCollection, int maxBranchCull)
{
	Timer timer("TrimSequences");
	printf("Trimming short branches: %u\n", maxBranchCull);
	unsigned numBranchesRemoved = 0;

	for (ISequenceCollection::iterator iter = seqCollection->begin();
			iter != seqCollection->end(); ++iter) {
		if (iter->second.deleted())
			continue;

		extDirection dir;
		// dir will be set to the trimming direction if the sequence can be trimmed
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
		// Sequence is trimmable, continue

		// This is a dead-end branch, check it for removal
		BranchRecord currBranch(dir, maxBranchCull);

		Kmer currSeq = iter->first;
		while(currBranch.isActive())
		{		
			// Get the extensions for this sequence, this function populates the extRecord structure
			ExtensionRecord extRec;
			int multiplicity = -1;
			bool success = seqCollection->getSeqData(currSeq, extRec, multiplicity);
			assert(success);
			(void)success;
			
			// process the extension record and extend the current branch, this function updates currSeq on successful extension
			processLinearExtensionForBranch(currBranch, currSeq, extRec, multiplicity);
		}
		
		// The branch has ended check it for removal, returns true if it was removed
		if(processTerminatedBranchTrim(seqCollection, currBranch))
		{
			numBranchesRemoved++;
		}
		seqCollection->pumpNetwork();
	}

	unsigned numSweeped = removeMarked(seqCollection);

	if (numBranchesRemoved > 0)
		PrintDebug(0, "Trimmed %u k-mer in %u branches\n",
				numSweeped, numBranchesRemoved);
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

// Process the extension for this branch for the trimming algorithm
// CurrSeq is the current sequence being inspected (the next member to be added to the branch). The extension record is the extensions of that sequence and
// multiplicity is the number of times that kmer appears in the data set
// After processing currSeq is unchanged if the branch is no longer active or else it is the generated extension
// If the parameter addKmer is true, add the k-mer to the branch.
bool processLinearExtensionForBranch(BranchRecord& branch,
		Kmer& currSeq, ExtensionRecord extensions, int multiplicity,
		bool addKmer)
{
	extDirection dir = branch.getDirection();
	extDirection oppDir = oppositeDirection(dir);
	
	if(branch.isTooLong())
	{
		// Check if the branch has extended past the max trim length.
		branch.terminate(BS_TOO_LONG);
		return false;
	} else if (extensions.dir[oppDir].isAmbiguous()) {
		// There is a reverse ambiguity to this branch, stop the branch without adding the current sequence to it
		branch.terminate(BS_AMBI_OPP);
		return false;
	}

	if (addKmer)
		branch.push_back(make_pair(currSeq,
					KmerData(multiplicity, extensions)));
	if (branch.isTooLong()) {
		branch.terminate(BS_TOO_LONG);
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
		PrintDebug(5, "Trimming %zu %s\n", branch.size(),
				branch.front().first.decode().c_str());
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
unsigned removeMarked(ISequenceCollection* pSC)
{
	Timer timer(__func__);
	unsigned count = 0;
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
		PrintDebug(1, "Removed %u marked k-mer\n", count);
	return count;
}

/** Assemble a contig.
 * @return the number of k-mer below the coverage threshold
 */
unsigned assembleContig(
		ISequenceCollection* seqCollection, FastaWriter* writer,
		BranchRecord& branch, unsigned id)
{
	assert(!branch.isActive());
	assert(branch.getState() == BS_NOEXT
			|| branch.getState() == BS_AMBI_SAME
			|| branch.getState() == BS_AMBI_OPP);

	// Assemble the contig.
	Sequence contig(branch);

	unsigned kmerCount = branch.calculateBranchMultiplicity();
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
unsigned assemble(ISequenceCollection* seqCollection,
		FastaWriter* fileWriter)
{
	Timer timer("Assemble");

	unsigned kmerCount = 0;
	unsigned contigID = 0;
	unsigned assembledKmer = 0;
	unsigned lowCoverageKmer = 0;
	unsigned lowCoverageContigs = 0;

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
			// singleton, output
			BranchRecord currBranch(SENSE, -1);
			currBranch.push_back(*iter);
			currBranch.terminate(BS_NOEXT);
			unsigned removed = assembleContig(seqCollection,
					fileWriter, currBranch, contigID++);
			assembledKmer += currBranch.size();
			if (removed > 0) {
				lowCoverageContigs++;
				lowCoverageKmer += removed;
			}
			continue;
		}
		assert(status == SC_ENDPOINT);

		BranchRecord currBranch(dir, -1);
		currBranch.push_back(*iter);
		Kmer currSeq = iter->first;
		extendBranch(currBranch, currSeq,
				iter->second.getExtension(dir));
		assert(currBranch.isActive());
		while(currBranch.isActive())
		{		
			// Get the extensions for this sequence, this function populates the extRecord structure
			ExtensionRecord extRec;
			int multiplicity = -1;
			bool success = seqCollection->getSeqData(currSeq, extRec, multiplicity);

			assert(success);
			(void)success;
			
			// process the extension record and extend the current branch, this function updates currSeq on successful extension
			processLinearExtensionForBranch(currBranch, currSeq, extRec, multiplicity);
		}
		
		if (currBranch.isCanonical()) {
			unsigned removed = assembleContig(seqCollection,
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
		printf("Found %u k-mer in %u contigs before removing "
				"low-coverage contigs\n", assembledKmer, contigID);
		printf("Removed %u k-mer in %u low-coverage contigs\n",
				lowCoverageKmer, lowCoverageContigs);
	} else {
		assert(assembledKmer <= kmerCount);
		unsigned circularKmer = kmerCount - assembledKmer;
		if (circularKmer > 0)
			printf("%u unassembled k-mer in circular contigs\n",
					circularKmer);
		printf("Assembled %u k-mer in %u contigs\n",
				assembledKmer, contigID);
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
		assert_open(histFile, opt::coverageHistPath);
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
