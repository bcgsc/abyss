#include "AssemblyAlgorithms.h"
#include "FastaReader.h"
#include "FastaWriter.h"
#include "Log.h"
#include "Options.h"
#include "ISequenceCollection.h"
#include "Timer.h"
#include "PackedSeqReader.h"
#include "PackedSeqWriter.h"
#include <cctype>
#include <cstdio>
#include <cstdlib>

namespace AssemblyAlgorithms
{

/** Return the kmer which are adjacent to this kmer. */
void generateSequencesFromExtension(const PackedSeq& currSeq, extDirection dir, SeqExt extension, PSequenceVector& outseqs)
{
	
	// Create the return structure
	PSequenceVector extensions;
	PackedSeq extSeq(currSeq);
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

/** Load packed sequences into the collection. */ 
static void loadPackedSequences(ISequenceCollection* seqCollection,
		std::string inFile)
{
	PackedSeqReader reader(inFile.c_str());

	unsigned count = 0;
	for (PSequenceVector seqs;
			reader.ReadSequences(seqs); seqs.clear()) {
		for (PSequenceVectorIterator iter = seqs.begin();
				iter != seqs.end(); iter++) {
			seqCollection->add(*iter);

			// Output the progress
			if (++count % 100000 == 0) {
				PrintDebug(1, "Read %u reads. ", count);
				seqCollection->printLoad();
			}
		}
		seqCollection->pumpNetwork();
	}
	PrintDebug(1, "Read %u reads. ", count);
	seqCollection->printLoad();

	if (count == 0)
		fputs("warning: input contains no sequences\n", stderr);
}

//
// Function to load sequences into the collection
//
void loadSequences(ISequenceCollection* seqCollection,
		std::string inFile)
{
	Timer timer("LoadSequences " + inFile);

	PrintDebug(0, "Reading `%s'\n", inFile.c_str());

	// Determine the input file type
	IFileReader* reader;
	if(inFile.find(".psq") != std::string::npos)
	{
		loadPackedSequences(seqCollection, inFile);
		return;
	}
	else if(inFile.find(".kmer") != std::string::npos) {
		seqCollection->load(inFile.c_str());
		return;
	}
	else
	{
		reader = new FastaReader(inFile.c_str());
	}

	unsigned count = 0, count_small = 0;
	for (SequenceVector seqs;
			reader->ReadSequences(seqs); seqs.clear()) {
		for (SequenceVectorIterator iter = seqs.begin();
				iter != seqs.end(); iter++) {
			int len = iter->length();
			if (opt::kmerSize > len) {
				count_small++;
				continue;
			}

			if (opt::rank <= 0
					&& count == 0 && seqCollection->count() == 0) {
				// Detect colour-space reads.
				seqCollection->setColourSpace(isdigit(iter->at(0)));
			} else {
				if (opt::colourSpace)
					assert(isdigit(iter->at(0)));
				else
					assert(isalpha(iter->at(0)));
			}

			for(int i = 0; i < len - opt::kmerSize  + 1; i++)
			{
				PackedSeq sub = iter->substr(i, opt::kmerSize);
				// Add the sequence to the sequence collection
				seqCollection->add(sub);
			}

			// Output the progress
			if (++count % 100000 == 0) {
				PrintDebug(1, "Read %u reads. ", count);
				seqCollection->printLoad();
			}
		}
		seqCollection->pumpNetwork();
	}
	PrintDebug(1, "Read %u reads. ", count);
	seqCollection->printLoad();

	unsigned count_nonacgt = reader->getNonACGT();
	delete reader;
	reader = 0;

	if (count_small > 0)
		fprintf(stderr, "warning: discarded %d sequences "
				"shorter than %d bases\n",
				count_small, opt::kmerSize);
	if (count_nonacgt > 0)
		fprintf(stderr, "warning: discarded %d sequences "
				"containing non-ACGT characters\n", count_nonacgt);

	if (count == 0)
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
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for (SequenceCollectionIterator iter = seqCollection->getStartIter();
			iter != endIter; ++iter) {
		if (iter->isFlagSet(SF_DELETE))
			continue;

		if (++count % 1000000 == 0)
			PrintDebug(1, "Generating adjacency: %d sequences\n", count);

		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			PackedSeq testSeq(*iter);
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

/** Remove all the extensions both from and to this sequence. */
static void removeExtensions(ISequenceCollection* seqCollection,
		const PackedSeq& seq, extDirection dir)
{
	removeExtensionsToSequence(seqCollection, seq, dir);
	seqCollection->clearExtensions(seq, dir);
}

/** Mark ambiguous branches and branches from palindromes for removal.
 * @return the number of branches marked
 */
unsigned markAmbiguous(ISequenceCollection* seqCollection)
{
	Timer timer(__func__);
	unsigned progress = 0;
	unsigned count = 0;
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		if(iter->isFlagSet(SF_DELETE))
			continue;

		if (++progress % 1000000 == 0)
			PrintDebug(1, "Splitting: %u sequences\n", progress);

		if (iter->isPalindrome()) {
			seqCollection->mark(*iter, SENSE);
			seqCollection->mark(*iter, ANTISENSE);
			count += 2;
			continue;
		}

		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			if (iter->getExtension(dir).isAmbiguous()
					|| iter->isPalindrome(dir)) {
				seqCollection->mark(*iter, dir);
				count++;
			}
		}
	}
	PrintDebug(0, "Marked %u ambiguous branches\n", count);
	return count;
}

/** Remove marked branches.
 * @return the number of branches removed
 */
unsigned splitAmbiguous(ISequenceCollection* pSC)
{
	Timer timer(__func__);
	unsigned count = 0;
	SequenceCollectionIterator end = pSC->getEndIter();
	for (SequenceCollectionIterator it = pSC->getStartIter();
			it != end; ++it) {
		if (pSC->checkFlag(*it, SF_DELETE))
			continue;
		for (extDirection sense = SENSE;
				sense <= ANTISENSE; ++sense) {
			if (pSC->isMarked(*it, sense)) {
				removeExtensions(pSC, *it, sense);
				count++;
			}
		}
		pSC->pumpNetwork();
	}
	PrintDebug(0, "Split %u ambiguous branches\n", count);
	return count;
}

int popBubbles(ISequenceCollection* seqCollection, int kmerSize)
{
	Timer timer("PopBubbles");
	int numPopped = 0;

	// Set the cutoffs
	const unsigned int expectedBubbleSize = 2*(kmerSize + 1);
	const unsigned int maxNumBranches = 3;
		
	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		// Skip sequences that have already been deleted	
		if(iter->isFlagSet(SF_DELETE))
		{		
			continue;
		}

		// Get the extensions for this sequence, this function populates the extRecord structure
		ExtensionRecord extRec;
		int multiplicity = -1;
		bool success = seqCollection->getSeqData(*iter, extRec, multiplicity);
		assert(success);
		(void)success;

		// Check for ambiguity
		for (extDirection dir = SENSE; dir <= ANTISENSE; ++dir) {
			if (extRec.dir[dir].isAmbiguous()) {
				// Found a potential bubble, examine each branch
				bool stop = false;
				
				// Create the branch group
				BranchGroup branchGroup(0, dir, maxNumBranches,
						*iter);
				initiateBranchGroup(branchGroup, *iter, extRec.dir[dir], multiplicity, expectedBubbleSize);

				// Disallow any further branching.
				unsigned numInitialBranches
					= branchGroup.getNumBranches();
				if (numInitialBranches <= maxNumBranches)
					branchGroup.setMaxNumBranches(numInitialBranches);
				else
					stop = true;
				
				// Iterate over the branches
				while(!stop)
				{
					size_t numBranches = branchGroup.getNumBranches();
					
					for(unsigned int j = 0; j < numBranches; ++j)
					{						
						// Get the extensions of this branch
						ExtensionRecord extRec;
						int multiplicity = -1;
						bool success = seqCollection->getSeqData(branchGroup.getBranch(j).getLastSeq(), extRec, multiplicity);
						assert(success);
						(void)success;
						
						processBranchGroupExtension(branchGroup, j, branchGroup.getBranch(j).getLastSeq(), extRec, multiplicity);
					}
					
					// At this point all branches should have the same length or one will be a noext
					branchGroup.updateStatus();
					
					// All branches have been extended one sequence, check the stop conditions
					BranchGroupStatus status = branchGroup.getStatus();
					
					// Check if a stop condition was met
					if(status == BGS_TOOLONG || status == BGS_LOOPFOUND || status == BGS_TOOMANYBRANCHES || status == BGS_NOEXT)
					{
						stop = true;
					}
					else if(status == BGS_JOINED)
					{
						static unsigned snpID;
						writeSNP(branchGroup, ++snpID);
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

	if (opt::snpFile != NULL)
		fflush(opt::snpFile);

	if (numPopped > 0)
		printf("Removed %d bubbles\n", numPopped);
	return numPopped;
}

//
// Populate a branch group with the inital branches from a sequence
//
void initiateBranchGroup(BranchGroup& group, const PackedSeq& seq, const SeqExt& extension, int /*multiplicity*/, size_t maxBubbleSize)
{
	
	// As the root sequence is not added to the branch, its multiplicity information is ignored.
	
	// Generate the vector of sequences that make up this branch
	PSequenceVector extSeqs;
	generateSequencesFromExtension(seq, group.getDirection(), extension, extSeqs);
	assert(extSeqs.size() > 1);
	uint64_t id = 0;
	
	for(PSequenceVector::iterator seqIter = extSeqs.begin(); seqIter != extSeqs.end(); ++seqIter)
	{
		// Create a new branch and add it to the group
		BranchRecord newBranch(group.getDirection(), maxBubbleSize);
		BranchRecord& addedBranch = group.addBranch(id, newBranch);
		
		// Add the sequence to the branch
		addedBranch.addSequence(*seqIter);
		
		id++;
	}	
}

//
// process an a branch group extension
//
bool processBranchGroupExtension(BranchGroup& group, size_t branchIndex, const PackedSeq& seq, ExtensionRecord extensions, int multiplicity)
{
	// Generate the extensions of the branch
	PSequenceVector branchExtSeqs;
	extDirection dir = group.getDirection();
	generateSequencesFromExtension(seq, dir, extensions.dir[dir], branchExtSeqs);
	
	// Set the multiplicity of the request sequence
	group.getBranch(branchIndex).setMultiplicity(seq, multiplicity);
	
	if(branchExtSeqs.size() == 1)
	{
		// single extension
		group.getBranch(branchIndex).addSequence(branchExtSeqs.front());
							
	}
	else if(branchExtSeqs.size() > 1)
	{
		// Start a new branch for the sequences [1..n]
		PSequenceVector::iterator seqIter = branchExtSeqs.begin() + 1;
		
		
		for(; seqIter != branchExtSeqs.end(); ++seqIter)
		{
			uint64_t newID = group.getNumBranches();
			
			// Start a new branch which is a duplicate of the current branch up to this point
			BranchRecord newBranch(group.getBranch(branchIndex));
			BranchRecord& addedBranch = group.addBranch(newID, newBranch);
			addedBranch.addSequence(*seqIter);
		}
		
		// Add the first sequence (index 0) to the current branch
		group.getBranch(branchIndex).addSequence(branchExtSeqs.front());
	}
	else
	{
	
		// this branch could not be extended, set a flag
		group.setNoExtension();
	}
	
	// Return whether the group is extendable
	return group.isExtendable();
}

/** Write SNP information to a file. */
void writeSNP(BranchGroup& group, unsigned id)
{
	FILE* fout = opt::snpFile;
	if (fout == NULL)
		return;

	unsigned selectedIndex = group.getBranchToKeep();
	char allele = 'A';

	BranchRecord& refBranch = group.getBranch(selectedIndex);
	Sequence refContig;
	refBranch.buildContig(refContig);
	fprintf(fout, ">%u%c %zu %u\n%s\n", id, allele++,
			refContig.length(),
			refBranch.getBranchMultiplicity(),
			refContig.c_str());

	unsigned numBranches = group.getNumBranches();
	for (unsigned i = 0; i < numBranches; ++i) {
		if (i == selectedIndex)
			continue;
		BranchRecord& currBranch = group.getBranch(i);
		Sequence contig;
		currBranch.buildContig(contig);
		fprintf(fout, ">%u%c %zu %u\n%s\n", id, allele++,
				contig.length(),
				currBranch.getBranchMultiplicity(),
				contig.c_str());
	}
}

//
// Collapse joined paths into a single path
//
void collapseJoinedBranches(ISequenceCollection* seqCollection, BranchGroup& group)
{
	assert(group.isAmbiguous(seqCollection));
	// a join was found, select a branch to keep and remove the rest
	size_t selectedIndex = group.getBranchToKeep();
	
	size_t numBranches = group.getNumBranches();
	BranchRecord& refRecord = group.getBranch(selectedIndex);
	PrintDebug(5, "Popping %zu %s\n", refRecord.getLength(),
				refRecord.getFirstSeq().decode().c_str());

	for(size_t i = 0; i < numBranches; ++i)
	{
		// Skip the branch that was selected to keep
		if(i == selectedIndex)
		{
			continue;
		}
		
		BranchRecord& currBranch = group.getBranch(i);
		BranchDataIter end = currBranch.getEndIter();
		for(BranchDataIter branchIter = currBranch.getStartIter();
				branchIter != end; ++branchIter)
		{
			/* As long as we're only popping simple bubbles, the
			 * sequence being removed cannot be in the reference
			 * sequence. By now, we've forgotten the multiplicity map
			 * used by BranchRecord::exists to save memory. */
			//assert(!refRecord.exists(*branchIter));
			removeSequenceAndExtensions(seqCollection, *branchIter);
		}
	}
	assert(!group.isAmbiguous(seqCollection));
}
//
// Remove a read and update the extension records of the sequences that extend to it
//
void removeSequenceAndExtensions(ISequenceCollection* seqCollection, const PackedSeq& seq)
{
	// This removes the reverse complement as well
	seqCollection->remove(seq);
	removeExtensionsToSequence(seqCollection, seq, SENSE);
	removeExtensionsToSequence(seqCollection, seq, ANTISENSE);
}

/** Remove all the extensions to this sequence. */
void removeExtensionsToSequence(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection dir)
{
	extDirection oppDir = oppositeDirection(dir);
	
	PackedSeq testSeq(seq);
	uint8_t extBase = testSeq.shift(dir);
	for(int i = 0; i < NUM_BASES; i++)
	{
		// don't bother checking if the extension exists, just remove it
		// if the sequence didnt have the extension, the operation will do nothing
		testSeq.setLastBase(dir, i);
		
		// remove the extension, this removes the reverse complement as well
		seqCollection->removeExtension(testSeq, oppDir, extBase);
	}	
}

/** The number of sequences that have been eroded. */
static unsigned g_numEroded;

/** Return the number of sequences that have been eroded. */
unsigned getNumEroded()
{
	unsigned numEroded = g_numEroded;
	g_numEroded = 0;
	PrintDebug(0, "Eroded %d tips\n", numEroded);
	return numEroded;
}

/** Consider the specified sequence for erosion.
 * @return the number of sequences eroded, zero or one
 */
unsigned erode(ISequenceCollection* c, const PackedSeq& seq)
{
	extDirection dir;
	SeqContiguity contiguity = checkSeqContiguity(c, seq, dir);
	if (contiguity == SC_INVALID || contiguity == SC_CONTIGUOUS)
		return 0;

	if (seq.getMultiplicity(SENSE) < 1
			|| seq.getMultiplicity(ANTISENSE) < 1) {
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
	Timer erodeEndsTimer("Erode sequences");
	assert(g_numEroded == 0);
	seqCollection->attach(erosionObserver);

	SequenceCollectionIterator endIter = seqCollection->getEndIter();
	for (SequenceCollectionIterator iter
			= seqCollection->getStartIter();
			iter != endIter; ++iter) {
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

//
// Check the adjacency of a sequence
//
SeqContiguity checkSeqContiguity(ISequenceCollection* seqCollection, const PackedSeq& seq, extDirection& outDir)
{
	if (seqCollection->checkFlag(seq, SF_DELETE))
		return SC_INVALID;
			
	bool child = seqCollection->hasChild(seq);
	bool parent = seqCollection->hasParent(seq);
	
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
	printf("Trimming short branches: %d\n", maxBranchCull);	
	unsigned numBranchesRemoved = 0;

	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		
		extDirection dir;
		// dir will be set to the trimming direction if the sequence can be trimmed
		SeqContiguity status = checkSeqContiguity(seqCollection, *iter, dir);

		if(status == SC_INVALID || status == SC_CONTIGUOUS)
		{
			continue;
		}
		else if(status == SC_ISLAND)
		{
			// remove this sequence, it has no extensions
			seqCollection->mark(*iter);
			continue;
		}
		// Sequence is trimmable, continue

		// This is a dead-end branch, check it for removal
		BranchRecord currBranch(dir, maxBranchCull);
					
		PackedSeq currSeq = *iter;
		
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
		PrintDebug(0, "Trimmed %u sequences in %u branches\n",
				numSweeped, numBranchesRemoved);
	return numBranchesRemoved;
}

//
// Process the extension for this branch for the trimming algorithm
// CurrSeq is the current sequence being inspected (the next member to be added to the branch). The extension record is the extensions of that sequence and
// multiplicity is the number of times that kmer appears in the data set
// After processing currSeq is unchanged if the branch is no longer active or else it is the generated extension
//
bool processLinearExtensionForBranch(BranchRecord& branch, PackedSeq& currSeq, ExtensionRecord extensions, int multiplicity)
{
	extDirection dir = branch.getDirection();
	extDirection oppDir = oppositeDirection(dir);
	
	if(branch.isTooLong())
	{
		// Check if the branch has extended past the max trim length.
		branch.terminate(BS_TOO_LONG);
	}
	else if(branch.hasLoop())
	{
		branch.terminate(BS_LOOP);
	} else if (extensions.dir[oppDir].isAmbiguous()) {
		// There is a reverse ambiguity to this branch, stop the branch without adding the current sequence to it
		branch.terminate(BS_AMBI_OPP);
	} else if (!extensions.dir[dir].hasExtension()) {
		// no extenstion, add the current sequence and terminate the branch
		branch.addSequence(currSeq, multiplicity);
		branch.terminate(branch.isTooLong() ? BS_TOO_LONG : BS_NOEXT);
	} else if (extensions.dir[dir].isAmbiguous()) {
		// this branch has an ambiguous extension, add the current sequence and terminate
		branch.addSequence(currSeq, multiplicity);
		branch.terminate(BS_AMBI_SAME);
	}
	else
	{
		// Add the sequence to the branch
		branch.addSequence(currSeq, multiplicity);
		
		// generate the new current sequence from the extension
		//printf("currseq: %s ", currSeq.decode().c_str());
		PSequenceVector newSeqs;

		generateSequencesFromExtension(currSeq, dir, extensions.dir[dir], newSeqs);
		assert(newSeqs.size() == 1);
		currSeq = newSeqs.front();
		//printf("newseq: %s \n", currSeq.decode().c_str());
	}
	
	return branch.isActive();
}

//
//
//
bool processTerminatedBranchTrim(ISequenceCollection* seqCollection, BranchRecord& branch)
{
	assert(!branch.isActive());
	if(branch.getLength() > 0 && branch.getState() != BS_TOO_LONG)
	{
		PrintDebug(5, "Trimming %zu %s\n", branch.getLength(),
					branch.getFirstSeq().decode().c_str());
		BranchDataIter endIter  = branch.getEndIter();
		for (BranchDataIter bIter = branch.getStartIter();
				bIter != endIter; bIter++)
			seqCollection->mark(*bIter);
		return true;
	}	
	else
	{
		return false;
	}
}

/** Remove all marked sequences.
 * @return the number of removed sequences
 */
unsigned removeMarked(ISequenceCollection* pSC)
{
	Timer timer(__func__);
	unsigned count = 0;
	SequenceCollectionIterator end = pSC->getEndIter();
	for (SequenceCollectionIterator it = pSC->getStartIter();
			it != end; ++it) {
		if (pSC->checkFlag(*it, SF_DELETE))
			continue;
		if (pSC->isMarked(*it)) {
			removeSequenceAndExtensions(pSC, *it);
			count++;
		}
		pSC->pumpNetwork();
	}
	if (count > 0)
		PrintDebug(1, "Removed %u marked sequences\n", count);
	return count;
}

/** Assemble a contig.
 * @return the number of k-mer below the coverage threshold
 */
unsigned assembleContig(
		ISequenceCollection* seqCollection, IFileWriter* writer,
		BranchRecord& branch, unsigned id)
{
	// Assemble the contig.
	Sequence contig;
	AssemblyAlgorithms::processTerminatedBranchAssemble(
			seqCollection, branch, contig);

	unsigned kmerCount = branch.calculateBranchMultiplicity();
	if (writer != NULL)
		writer->WriteSequence(contig, id, kmerCount);

	// Remove low-coverage contigs.
	float coverage = (float)kmerCount / branch.getLength();
	BranchDataIter end = branch.getEndIter();
	if (opt::coverage > 0 && coverage < opt::coverage) {
		for (BranchDataIter it
				= branch.getStartIter();
				it != end; ++it)
			seqCollection->remove(*it);
		return branch.getLength();
	}
	return 0;
}

/** Assemble contigs.
 * @return the number of contigs assembled
 */
unsigned assemble(ISequenceCollection* seqCollection,
		IFileWriter* fileWriter)
{
	Timer timer("Assemble");

	unsigned contigID = 0;
	unsigned lowCoverageKmer = 0;
	unsigned lowCoverageContigs = 0;

	SequenceCollectionIterator endIter  = seqCollection->getEndIter();
	for(SequenceCollectionIterator iter = seqCollection->getStartIter(); iter != endIter; ++iter)
	{
		
		extDirection dir;
		// dir will be set to the trimming direction if the sequence can be trimmed
		SeqContiguity status = checkSeqContiguity(seqCollection, *iter, dir);

		if(status == SC_INVALID || status == SC_CONTIGUOUS)
		{
			continue;
		}
		else if(status == SC_ISLAND)
		{
			// singleton, output
			BranchRecord currBranch(SENSE, -1);
			currBranch.addSequence(*iter, iter->getMultiplicity());
			currBranch.terminate(BS_NOEXT);
			unsigned removed = assembleContig(seqCollection,
					fileWriter, currBranch, contigID++);
			if (removed > 0) {
				lowCoverageContigs++;
				lowCoverageKmer += removed;
			}
			continue;
		}

		// The sequence is an endpoint, begin extending it
		// Passing -1 into the branch will disable the length check
		BranchRecord currBranch(dir, -1);
					
		PackedSeq currSeq = *iter;
		
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
			if (removed > 0) {
				lowCoverageContigs++;
				lowCoverageKmer += removed;
			}
		}

		seqCollection->pumpNetwork();
	}

	if (opt::coverage > 0) {
		printf("Found %u contigs before removing "
				"low-coverage contigs\n", contigID);
		printf("Removed %u k-mer in %u low-coverage contigs\n",
				lowCoverageKmer, lowCoverageContigs);
	} else
		printf("Assembled %u contigs\n", contigID);
	return contigID;
}

//
//
//
void processTerminatedBranchAssemble(
		ISequenceCollection* /*seqCollection*/,
		const BranchRecord& branch, Sequence& outseq)
{
	assert(!branch.isActive());
	//printf("	branch has size: %d\n", branchElements.size());
	
	// the only acceptable condition for the termination of an assembly is a noext or a loop
	assert(branch.getState() == BS_NOEXT || branch.getState() == BS_LOOP);
	
	// Assemble the contig
	branch.buildContig(outseq);
}

// Write the sequences out to a file
//
void outputSequences(const char* filename, ISequenceCollection* pSS)
{
	FastaWriter writer(filename);
	SequenceCollectionIterator endIter  = pSS->getEndIter();
	int64_t count = 0;
	for(SequenceCollectionIterator iter = pSS->getStartIter(); iter != endIter; ++iter)
	{
		if(!pSS->checkFlag(*iter, SF_DELETE))
		{
			writer.WriteSequence(iter->decode(), count,
					iter->getMultiplicity());
			count++;
		}
	}	
}

// Write the sequences out to a file
//
void outputPackedSequences(const char* filename, ISequenceCollection* pSS)
{
	PackedSeqWriter writer(filename);
	SequenceCollectionIterator endIter  = pSS->getEndIter();
	int64_t count = 0;
	for(SequenceCollectionIterator iter = pSS->getStartIter(); iter != endIter; ++iter)
	{
		
		if(!pSS->checkFlag(*iter, SF_DELETE))
		{
			writer.WriteSequence(*iter);
			count++;
		}
	}	
}

};
