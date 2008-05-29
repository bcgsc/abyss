#include <stdio.h>
#include <math.h>
#include <iostream>
#include "Scaffold.h"
#include "Sequence.h"
#include "FastaWriter.h"
#include "SequenceCollectionHash.h"
#include "Aligner.h"
#include "FastaReader.h"
#include "ParentTree.h"
#include "DirectedGraph.h"
#include "Options.h"
#include "ContigData.h"

using namespace std;

int main(int argc, char** argv)
{	
	if(argc < 4)
	{
		printf("usage: Scaffold <final reads file> <readsFile> <contigFile> <readLen> <assembly kmer size> [alignment file]\n");
		exit(1);
	}
	
	std::string finalReadsFile(argv[1]);
	std::string readsFile(argv[2]);
	std::string contigFile(argv[3]);
	int readLen = atoi(argv[4]);
	int kmer = atoi(argv[5]);
	
	// kinda hacky
	opt::kmerSize = kmer;
	opt::readLen = readLen;
	
	// Should we read in the alignment hints file?
	//bool readAlignments = false;
	//std::string alignmentsFile = "";
	
	std::string pairsCacheFile = "";
	std::string alignmentCacheFile = "";
	
	if(argc > 6)
	{
		pairsCacheFile = argv[6];
	}
	
	if(argc > 7)
	{
		alignmentCacheFile = argv[7];
	}
	
	Scaffold scaffold(finalReadsFile, readsFile, contigFile, readLen, kmer, pairsCacheFile, alignmentCacheFile);
}

//
//
//
Scaffold::Scaffold(std::string finalReadsFile, std::string readsFile, std::string contigFile, int readLen, int kmer, std::string pairsCacheFile, std::string /*alignCacheFile*/) : m_readLen(readLen), m_kmer(kmer)
{
	// Read in the pairs
	//ReadAlignments(alignFile);
	
	// Read in the adjacency graph
	LoadAdjacency(finalReadsFile);
			
	// Read in the contigs
	PairedAlgorithms::ReadContigs(contigFile, m_contigMap);

	// Read in the sequencing reads
	if(pairsCacheFile.empty())
	{
		ReadSequenceReads(readsFile);	
		LoadPairsRecord(m_readVec, m_kmer);
		m_pairRec.serialize("PairsCache.bin");		
	}
	else
	{
		m_pairRec.unserialize(pairsCacheFile);	
	}
	
	/*
	// read in the alignments
	if(alignCacheFile.empty())
	{
		printf("Generating alignments from reads file\n");
		GenerateAlignmentCache(m_pSC, m_contigMap, m_alignCache);
		m_alignCache.serialize("alignCache.bin");
		
	}
	else
	{
		printf("Reading the alignments from file %s\n", alignCacheFile.c_str());
		m_alignCache.unserialize(alignCacheFile);
	}
	return;
	*/
	merge4();
	exit(1);
	

	// Generate the initial alignment database
	
	

	

	
	printf("Done reading\n");
	// Generate the empirical distribution
	
	GenerateStatistics();

	merge2();
	printf("Done all merges, outputting\n");
	/*
	bool stop = false;
	while(!stop)
	{
		bool anyMerge = false;
		for(CMIter iter = m_contigMap.begin(); iter != m_contigMap.end(); iter++)
		{
			if(!iter->second.merged && iter->second.seq.length() > 200)
			{
				if(AttemptMerge(iter->first))
				{
					anyMerge = true;
				}
			}
		}
		
		if(!anyMerge)
		{
			stop = true;
		}
	}
	*/
	// Output the merged contigs
	ofstream mergeFile("merged.fa");
	for(CMIter iter = m_contigMap.begin(); iter != m_contigMap.end(); iter++)
	{
		if(!iter->second.merged)
		{
			mergeFile << ">" << iter->first << " " << iter->second.seq.length() << endl << iter->second.seq << endl; 
		}
	}
	mergeFile.close();
}

Scaffold::~Scaffold()
{
	delete m_pSC;	
}



void Scaffold::merge2()
{
	ofstream mergeFile("merged2.fa");

	// Examine each end of the contig to see if it can be extended by pairs
	for(CMIter iter = m_contigMap.begin(); iter != m_contigMap.end(); iter++)
	{
		std::set<ContigID> mergeIDs;
		if(!iter->second.merged && iter->second.seq.length() > 200)
		{
			for(int toggle = 0; toggle <= 1; toggle++)
			{
				bool stop = false;
				
				if(toggle == 1)
				{
					printf("trying rev comp\n");
					iter->second.seq = reverseComplement(iter->second.seq);
				}
				
				while(!stop)
				{
					printf("trying %s for merge\n", iter->first.c_str());
					
					// Get the last k bases from the end
					Sequence& seq = iter->second.seq;
					
					int start = seq.length() - m_kmer;
					PackedSeq endMer(seq.substr(start, m_kmer));
					
					// Build the tree from this endmer
					int lookForwardDepth = 200;
					int lookBackDepth = 200;
	
					ParentTree parentTree(endMer, SENSE, m_pSC, lookForwardDepth);
	
					// Get pairs
					int numNodesInBranch = seq.length();
					size_t firstIndex;
					size_t lastIndex;
					
					lastIndex = numNodesInBranch - m_kmer;
									
									
					if((int)lastIndex - lookBackDepth < 0)
					{
						firstIndex = 0;
					}
					else
					{
						firstIndex = lastIndex - lookBackDepth;
					}
													
					PSeqSet children = parentTree.getRootChildren();					
					const double minSupport = 10.0f;
					
					PSequenceVector allPairs;
					
					int numGood = 0;
					
					for(size_t idx = firstIndex; idx < lastIndex; ++idx)
					{
						assert(idx + m_kmer < seq.length());
						PackedSeq firstHalf(seq.substr(idx, m_kmer));
						
						double weight = 1.0f;
						PSequenceVector seqPairs = m_pairRec.getPairs(firstHalf);
						
						
						for(PSequenceVector::iterator pairIter = seqPairs.begin(); pairIter != seqPairs.end(); pairIter++)
						{
							// Check if this pair is supporting multiple branches
							// If it is, do not add the pairs to the set
							parentTree.addScore(*pairIter, weight);
						}
						
						
						// Check the number of branches supported by the pairs of this seq
						int numSupported = 0;
						int suppIdx = 0;
						int currIdx = 0;
						for(PSeqSet::iterator testIter = children.begin(); testIter != children.end(); ++testIter)
						{
							double currScore = parentTree.getScore(*testIter).subtreeScore;	
							if(currScore > 0)
							{
								suppIdx = currIdx;
								numSupported++;
							}
							currIdx++;
						}					
	
						if(numSupported == 1)
						{
							for(PSequenceVector::iterator pairIter = seqPairs.begin(); pairIter != seqPairs.end(); pairIter++)
							{
								allPairs.push_back(*pairIter);
							}
							numGood++;
						}
						parentTree.resetScore();
					}
					
					printf("%d of %zd sequences are usable\n",
							numGood, lastIndex - firstIndex);
					
					int count = 0;
					for(PSequenceVector::iterator pairIter = allPairs.begin(); pairIter != allPairs.end(); pairIter++)
					{
						count++;
						parentTree.addScore(*pairIter, 1.0f);
					}
					
					parentTree.print(10);
					
					parentTree.printRootChildren();
			
					PSeqSet::iterator bestIter;
					
					int numSupported = 0;
					for(PSeqSet::iterator testIter = children.begin(); testIter != children.end(); ++testIter)
					{
						double currScore = parentTree.getScore(*testIter).subtreeScore;	
						if(currScore > minSupport)
						{
							bestIter = testIter;
							numSupported++;
						}
					}
					
					if(numSupported == 1)
					{
						//perform a merge
						
						// find the contig with the specified kmer as a start
						ContigID slaveID;
						Sequence found = findContig(*bestIter, slaveID);
						if(!found.empty() && slaveID != iter->first && mergeIDs.find(slaveID) == mergeIDs.end())
						{
							printf("merging with %s\n", found.c_str());
							
							// perform merge
							seq.append(found.substr(m_kmer - 1));
							
							// replace
							// Replace the contig with the merged product
							//m_contigMap[contigID].seq = merged;
							
							// Erase the existance of the slave contig
							//m_contigMap.erase(slaveID);
							iter->second.super = true;	
							m_contigMap[slaveID].merged = true;
							mergeIDs.insert(slaveID);				
						}
						else
						{
							// only append if the seq doesnt exist 
							if(seq.find(bestIter->decode()) == string::npos)
							{
								seq.append(1, bestIter->getLastBase());
							}
							else
							{
								stop = true;
							}
							printf("appending last seq only\n");	
						}
						
						printf("merged %s\n", seq.c_str());
						
					}
					else
					{
						stop = true;
					}
				}
			}
		}
		
		if(iter->second.super)
		{
			mergeFile << ">" << iter->first << " " << iter->second.seq.length() << endl << iter->second.seq << endl; 
		}
	}
	mergeFile.close();
}

struct CompareBranches
{
	bool operator()(const BranchRecord& r1, const BranchRecord& r2) { return r1.getLength() < r2.getLength(); }
};

void Scaffold::merge3()
{
	Timer timer("MergeAssemble");
	FastaWriter* fileWriter = new FastaWriter("merge3contigs.fa");
	FastaWriter* compWriter = new FastaWriter("componentContigs.fa");
	
	// generate the master list of starting sequences
	contigStartVec starts;
	PairedAlgorithms::generateStartList(m_pSC, starts);
	
	printf("discovered %zu starts\n", starts.size());	
	int numAssembled = 0;
	int contigID = 0;

	std::priority_queue<BranchRecord, std::vector<BranchRecord>, CompareBranches> branchPQ;
	
	for(contigStartVec::iterator iter = starts.begin(); iter != starts.end(); iter++)
	{	
		if(m_pSC->checkFlag(iter->seq, SF_SEEN))
		{
			continue;
		}
		
		printf("assembing contig %d\n", contigID);
		
		PackedSeq currSeq = iter->seq;		
		if(iter->state == CSS_ISLAND)
		{
			// singleton, output
			BranchRecord currBranch(SENSE, -1);
			currBranch.addSequence(currSeq);
			currBranch.terminate(BS_NOEXT);
			
			// this will set the seen flag
			Sequence contig;
			AssemblyAlgorithms::processTerminatedBranchAssemble(m_pSC, currBranch, contig);
			compWriter->WriteSequence(contig, contigID, 0);
			contigID++;
			branchPQ.push(currBranch);
			continue;
		}
		
		extDirection dir = iter->dir;

		// The sequence is an endpoint, begin extending it
		// Passing -1 into the branch will disable the length check
		BranchRecord currBranch(dir, -1);
		
		while(currBranch.isActive())
		{		
			// Get the extensions for this sequence, this function populates the extRecord structure
			ExtensionRecord extRec;
			int multiplicity = -1;
			bool success = m_pSC->getSeqData(currSeq, extRec, multiplicity);

			assert(success);
			(void)success;
			
			
			// process the extension record and extend the current branch, this function updates currSeq on successful extension
			PairedAlgorithms::processNonlinearExtensionForBranch(m_pSC, &m_pairRec, currBranch, currSeq, extRec, multiplicity);
		}
		
		// hack
		currBranch.terminate(BS_NOEXT);

		// this will set the seen flag
		Sequence contig;
		AssemblyAlgorithms::processTerminatedBranchAssemble(m_pSC, currBranch, contig);
		compWriter->WriteSequence(contig, contigID, 0);
		contigID++;	
		branchPQ.push(currBranch);			
		numAssembled++;

		m_pSC->pumpNetwork();
	}
	
	// wipe the seen flag
	m_pSC->wipeFlag(SF_SEEN);
	
	// output the branches in terms of priority
	while(!branchPQ.empty())
	{
		const BranchRecord& curr = branchPQ.top();
		printf("branch len: %zu\n", curr.getLength());
		
		bool unique = false;
		
		// Output the top branch if all its sequences are not seen (or should this be if ANY sequence is not seen?)
		for(size_t idx = 0; idx < curr.getLength(); ++idx)
		{
			if(!m_pSC->checkFlag(curr.getSeqByIndex(idx), SF_SEEN))
			{
				unique = true;
				break;
			}
		}
		
		if(unique)
		{
			// output the branch
			Sequence contig;
			AssemblyAlgorithms::processTerminatedBranchAssemble(m_pSC, curr, contig);
			
			// Output the contig
			fileWriter->WriteSequence(contig, contigID, 0);
			contigID++;	
		}
		
		branchPQ.pop();
	}
	
	delete fileWriter;
	delete compWriter;
	
	printf("num assembled: %d\n", numAssembled);	
}

void Scaffold::merge4()
{
	printf("Building graph\n");
	
	// Create a graph of all the contigs
	DirectedGraph<ContigID, ContigData> contigGraph;
	
	// Add all the vertices
	for(ContigMap::iterator contigIter = m_contigMap.begin(); contigIter != m_contigMap.end(); ++contigIter)
	{
		ContigData data(contigIter->second.seq, m_kmer);
		contigGraph.addVertex(contigIter->first, data);
	}
	
	// Generate a k-mer -> contig lookup table for all the contig ends
	std::map<PackedSeq, ContigID> contigLUT;
	for(ContigMap::iterator contigIter = m_contigMap.begin(); contigIter != m_contigMap.end(); ++contigIter)
	{
		Sequence& contigSequence = contigIter->second.seq;
		const unsigned numEnds = 2;
		PackedSeq seqs[numEnds];
		seqs[0] = PackedSeq(contigSequence.substr(contigSequence.length() - m_kmer, m_kmer)); //SENSE
		seqs[1] = PackedSeq(contigSequence.substr(0, m_kmer)); // ANTISENSE
		
		size_t numToAdd = (seqs[0] != seqs[1]) ? 2 : 1;
		
		for(unsigned idx = 0; idx < numToAdd; idx++)
		{	
			// insert sequences into the table
			contigLUT[seqs[idx]] = contigIter->first;
		}
	}
	
	// Build the edges
	for(ContigMap::iterator contigIter = m_contigMap.begin(); contigIter != m_contigMap.end(); ++contigIter)
	{
		// Generate edges to/from this node
		
		// Since two contigs are not necessarily built from the same strand, two contigs can both have OUT nodes pointing to each other
		// this situation will get cleaned up when the links are resolved/merged
		Sequence& contigSequence = contigIter->second.seq;
		const unsigned numEnds = 2;
		PackedSeq seqs[numEnds];
		seqs[0] = PackedSeq(contigSequence.substr(contigSequence.length() - m_kmer, m_kmer)); //SENSE
		seqs[1] = PackedSeq(contigSequence.substr(0, m_kmer)); // ANTISENSE

		ExtensionRecord extRec;
		int multiplicity;
				
		for(unsigned idx = 0; idx < numEnds; idx++)
		{
			PackedSeq& currSeq = seqs[idx];
			extDirection dir;
			dir = (idx == 0) ? SENSE : ANTISENSE;
						
			m_pSC->getSeqData(currSeq, extRec, multiplicity);
			
			// Generate the links
			PSequenceVector extensions;
			AssemblyAlgorithms::generateSequencesFromExtension(currSeq, dir, extRec.dir[dir], extensions);
			
			for(PSequenceVector::iterator iter = extensions.begin(); iter != extensions.end(); ++iter)
			{
				// Get the contig this sequence maps to
				bool foundEdge = false;
				for(size_t compIdx = 0; compIdx <= 1; ++compIdx)
				{
					bool reverse = (compIdx == 1);
					PackedSeq testSeq;
					if(reverse)
					{
						testSeq = reverseComplement(*iter);
					}
					else
					{
						testSeq = *iter;
					}
				
					std::map<PackedSeq, ContigID>::iterator cLUTIter;
					cLUTIter = contigLUT.find(testSeq);
					if(cLUTIter != contigLUT.end())
					{						
						contigGraph.addEdge(contigIter->first, cLUTIter->second, dir, reverse);
						foundEdge = true;
					}
				}
				
				// it should ALWAYS be found since all sequences in the data set must belong to a contig		
				assert(foundEdge);
			}
		}

	}
	
	size_t numVert = contigGraph.getNumVertices();
	size_t numEdges = contigGraph.countEdges(); // SLOW
	printf("BEFORE: num vert: %zu num edges: %zu\n", numVert, numEdges);
		
	// Generate the initial alignment DB
	AlignmentCache alignDB;
	
	// Create the contig data function object
	ContigDataFunctions dataFunctor(m_kmer, &alignDB);	
	
	printf("Pre-validating graph\n");	
	contigGraph.validate(dataFunctor);

	FastaWriter* pWriter;
	pWriter = new FastaWriter("BeforeContigs.fa");
	
	contigGraph.iterativeVisit(ContigDataOutputter(pWriter));

	delete pWriter;
	pWriter = 0;

	// Populate the database
	{
		Timer t("GenDB");
		contigGraph.iterativeVisit(DBGenerator(&alignDB));
	}	
	
	bool stop = false;
	while(!stop)
	{
		size_t numMerged = contigGraph.removeTransitivity(dataFunctor);
		printf("trans reduce performed %zu merges\n", numMerged);
			
		if(numMerged == 0)
		{
			stop = true;
		}
		else
		{
			numVert = contigGraph.getNumVertices();
			numEdges = contigGraph.countEdges(); // SLOW
			printf("NOW: num vert: %zu num edges: %zu\n", numVert, numEdges);					
		}
	}
	
	AlignmentCache db2;
	contigGraph.iterativeVisit(DBGenerator(&db2));	
	
	
	printf("Post validating graph\n");
	contigGraph.validate(dataFunctor);

	size_t maxComponentLength = 500;
	PairedResolver resolver(m_kmer, maxComponentLength, &m_pairRec, &alignDB);
	
	SequenceDataCost dataCost;
	
	//contigGraph.reducePaired(maxComponentLength, dataCost, resolver, dataMerger);

	contigGraph.attemptResolve("40", SENSE, maxComponentLength, dataCost, resolver, dataFunctor);
	
	contigGraph.printVertex("40", true);
	
	return;
	
	numVert = contigGraph.getNumVertices();
	numEdges = contigGraph.countEdges(); // SLOW
	printf("AFTER: num vert: %zu num edges: %zu\n", numVert, numEdges);	

	pWriter = new FastaWriter("AfterSimplificationContigs.fa");
	
	contigGraph.iterativeVisit(ContigDataOutputter(pWriter));
	delete pWriter;
	pWriter = 0;
		
	printf("Attemping to resolve using pairs\n");
	

	
	pWriter = new FastaWriter("ResolvedContigs.fa");
	contigGraph.iterativeVisit(ContigDataOutputter(pWriter));
	delete pWriter;	
	contigGraph.validate(dataFunctor);
}

Sequence Scaffold::findContig(PackedSeq start, ContigID& id)
{
	for(CMIter iter = m_contigMap.begin(); iter != m_contigMap.end(); iter++)
	{
		if(iter->second.super)
		{
			continue;
		}
		
		Sequence test = iter->second.seq;
		
		// Get the first k bases
		PackedSeq sub1 = test.substr(0, m_kmer);
		
		// And the last k bases
		PackedSeq sub2 = test.substr(test.length() - m_kmer, m_kmer);
		
		if(start == sub1)
		{
			id = iter->first;
			return test;
		}
		else if(reverseComplement(start) == sub2)
		{
			id = iter->first;
			return reverseComplement(test);
		}
	}
	
	Sequence empty;
	return empty;
}

//
// Generate the initial alignments for the sequences
//
void Scaffold::GenerateAlignments(PSequenceVector& /*seqs*/, ContigMap& /*contigs*/)
{
	assert(false);
	/*
	printf("Creating DB\n");	
	Aligner aligner(m_kmer);
	aligner.CreateDatabase(m_contigMap);
	
	int id = 0;
	for(PSequenceVectorIterator iter = seqs.begin(); iter != seqs.end(); ++iter)
	{	
		// get the alignments
		AlignmentVector ar = aligner.GetAlignments(*iter);
		
		for(AlignmentVector::iterator arIter = ar.begin(); arIter != ar.end(); arIter++)
		{
			ReadAlign ra = BuildReadAlign(id, arIter->contig, arIter->start, arIter->isRC);
			
			// Set up all the mappings
			m_alignMap[id].push_back(ra);
			
			m_contigReadMap[arIter->contig].insert(id);
		}
	
		id++;
	}
	*/
}

//
// Get all the linkages between this contig and its pairs
//
bool Scaffold::AttemptMerge(ContigID contigID)
{
	bool mergedOccured = false;
	ContigPairVecMap cpvMap;
	// Get all the pair alignments for this contig
	// It will generate a mapping of contig pairs (c1,c2) to a vector of the pairs supporting the contigs
	// This function populates the cpvMap structure
	GenerateUniquePairAlignments(contigID, cpvMap);
	
	cout << "Contig " << contigID << " is linked to: " << endl;
	
	// Find the linkages of this contig on the left and right
	LinkVec linkages[2];
	
	for(CPVMIter cpvmIter = cpvMap.begin(); cpvmIter != cpvMap.end(); cpvmIter++)
	{		
		if(!m_contigMap[cpvmIter->first].repetitive)
		{
			ContigLinkage link = GenerateLinkage(contigID, cpvmIter->first, cpvmIter->second);
			
			if(!link.noLink)
			{
				printf("	link[%d]: %s at %d\n", link.order, link.slaveID.c_str(), link.distance);
				linkages[link.order].push_back(link);
			}
		}
	}
	
	for(int i = 0; i <= 1; i++)
	{
		printf("Contig %s (%zu bp) has %zu linkages on the %d side\n", contigID.c_str(), m_contigMap[contigID].seq.length(), linkages[i].size(), i);

		// Choose the best link
		bool hasBestLink = false;
		LinkIter bestLink;
		int bestDistance = -m_kmer;
		
		const int REQUIRED_PAIRS = 15;
		
		for(LinkIter iter = linkages[i].begin(); iter != linkages[i].end(); iter++)
		{
			printf("	link %s has %d pairs [%d %zu]\n", iter->slaveID.c_str(), iter->numPairs, iter->distance, m_contigMap[iter->slaveID].seq.length());
			if(iter->numPairs > REQUIRED_PAIRS)
			{
				if(iter->distance > bestDistance)
				{
					bestLink = iter;
					bestDistance = iter->distance;
					hasBestLink = true;
				}
			}
		}		
		
		// Is there a good link for this contig?
		if(hasBestLink)
		{
			printf("Best link is %s at %d bp away(pairs %d)\n", bestLink->slaveID.c_str(), bestLink->distance, bestLink->numPairs);
			
			bool consistent = CheckConsistency(*bestLink, linkages[i]);
			printf("Best link is %s\n", consistent ? "consistent" : "inconsistent");
			
			if(consistent)
			{
	
				// Get the sequences
				Sequence contig0 = m_contigMap[bestLink->masterID].seq;
				Sequence contig1 = m_contigMap[bestLink->slaveID].seq;
	
				// Create pointers to the sequences in the correct order
				Sequence* leftContig;
				Sequence* rightContig;
				
				// get the pairs of the sequence that can potentially fill the gap
				bool contig0Comp;
				bool contig1Comp;
				
				if(bestLink->order == CORDER_LEFT)
				{
					contig0Comp = false;
					contig1Comp = true;
					
					// set the pointers
					leftContig = &contig0;
					rightContig = &contig1;
				}
				else
				{
					contig0Comp = true;
					contig1Comp = false;		
					
					// set the pointers
					leftContig = &contig1;
					rightContig = &contig0;				
				}
				
				// If contig1 is in the RC direction, reverse the direction of the required pairs
				if(bestLink->orientation == CORIEN_OPP)
				{
					contig1Comp = !contig1Comp;
					contig1 = reverseComplement(contig1);
				}
	
				// Reverse assemble based on adjacency info		
				
				// Build the tree fanning out from the last kmer of the master contig
				int lookForwardDepth = bestDistance + 40;
				PackedSeq leftSeq(leftContig->substr(leftContig->length() - (m_kmer), m_kmer));
				PackedSeq rightSeq(rightContig->substr(0, m_kmer));
				
				printf("Assembling from %s to %s\n", leftSeq.decode().c_str(), rightSeq.decode().c_str());
				
				ParentTree parentTree(leftSeq, SENSE, m_pSC, lookForwardDepth);
				
				std::vector<PSequenceVector> outpaths;
				
				// empty to start
				PSequenceVector inseq;
				
				size_t max_paths = 20;
				parentTree.buildBackwards(rightSeq, leftSeq, inseq, outpaths, lookForwardDepth, max_paths);
				printf("Found %zu paths\n", outpaths.size());
				Sequence merged;
				
				// Choose best
				if(outpaths.size() > 0 && outpaths.size() <= max_paths)
				{				
					for(std::vector<PSequenceVector>::iterator pathIter = outpaths.begin(); pathIter != outpaths.end(); pathIter++)
					{
						// reverse the vector
						std::reverse(pathIter->begin(), pathIter->end());
						printf("Received %zu seqs\n", pathIter->size());
						// shatter the right contig and append it
						for(size_t i = 0; i < rightContig->length() - m_kmer; i++)
						{
							pathIter->push_back(PackedSeq(rightContig->substr(i, m_kmer)));
						}
						
						Sequence tempMerged = *leftContig;
						// Build the insert
						for(PSequenceVector::iterator insIter = (pathIter->begin() + 1); insIter != pathIter->end(); insIter++)
						{
							tempMerged.append(1, insIter->getLastBase());
						}
						
						printf("Merged contig: %s\n", tempMerged.c_str());
						merged = tempMerged;
					}
					
					// Perform the merge and update the distance
					int c0Len = m_contigMap[contigID].seq.length();
					int c1Len = m_contigMap[bestLink->slaveID].seq.length();
							
					// Update the reads for the master contig
					int offset = 0;
					if(bestLink->order == CORDER_RIGHT)
					{
						// Offset the master contig pairs by the increase in distance
						offset = merged.length() - c0Len;
					}
					
	
					UpdateMasterReads(contigID, offset, m_contigMap[contigID].seq, merged);
							
					// Update the slave record
					offset = 0;
					if(bestLink->order == CORDER_LEFT)
					{
						offset = merged.length() - c1Len;
					}
					
					bool isFlipped = bestLink->orientation == CORIEN_OPP;
					UpdateSlaveReads(bestLink->slaveID, contigID, offset, isFlipped, m_contigMap[bestLink->slaveID].seq, merged);
					
					// Replace the contig with the merged product
					m_contigMap[contigID].seq = merged;
					
					// Erase the existance of the slave contig
					m_contigMap.erase(bestLink->slaveID);
					
					mergedOccured = true;
					
					// Realign the pairs of the contig's reads
					RealignContigPairs(contigID);
					
					/*
					// Mark the intermediate links as repetitive so they are no longer considered
					std::set<ContigID> exclusion;
					for(LinkIter iter = linkages[i].begin(); iter != linkages[i].end(); iter++)
					{
						if(iter->slaveID != bestLink->slaveID)
						{
							m_contigMap[iter->slaveID].repetitive = true;
						}
					}
					*/
					
					printf("MERGED (%zu): %s\n", m_contigMap[contigID].seq.length(), m_contigMap[contigID].seq.c_str());
				}
			}
			
			
		}
	}

	return mergedOccured;
}

//
// Merge two contigs using the information provided in the link data structure
//
int Scaffold::Merge(Sequence& leftContig, Sequence& rightContig, int distance, Sequence& merged)
{
	printf("Merging at: %d\n", distance);	
	// Check if the contigs (potentially overlap)
	if(distance < 0)
	{
		distance = abs(distance);
		// Generate the merged sequence
		merged = leftContig;
		
		// Add the substring of the second sequence
		merged.append(rightContig.substr(distance, rightContig.length() - distance));

		//printf("ALIGN DIST: %d\n", alignDistance);
	}
	else
	{
		// Positive distance, fill the gap with "N"
		merged = leftContig;
		merged.append(distance, 'N');
		merged.append(rightContig);
	}
	return distance;
}

//
// Align contigs using the input position as a guess to the alignment
// Precondition: Contigs are from the same strand
// The estimated position is relative to the beginning of the left contig 
int Scaffold::alignContigs(const Sequence& leftContig, const Sequence& rightContig, int estimate, int range, int& retScore)
{
	int bestScore = -1000;
	int start = estimate - range;
	
	if(start < 0)
	{
		start = 0;
	}
	
	int stop = estimate + range;
	int leftLength = leftContig.length();
	int rightLength = rightContig.length();
	
	if(stop > leftLength)
	{
		stop = leftContig.length();
	}
	
	int bestI = start;	
	printf("Aligning: %s with %s estimated to be %d\n", leftContig.c_str(), rightContig.c_str(), estimate);
	
	for(int i = start; i < stop; i++)
	{
		
		// Compute amount of overlap between the contigs
		int overlap = leftContig.length() - i;
		
		// if the amount of overlap is longer than the right contig, trim the amount of overlap (the contig is a strict subset)
		if(overlap > rightLength)
		{
			overlap = rightContig.length();
		}
		
		// Get the subsequences
		Sequence leftsub = leftContig.substr(i, overlap);
		Sequence rightsub = rightContig.substr(0, overlap);
		
		assert(leftsub.length() == rightsub.length());
		
		int matched = 0;
		int unmatched = 0;
		for(int j = 0; j < overlap; j++)
		{
			if(leftsub[j] == rightsub[j])
			{
				matched++;
			}
			else
			{
				unmatched++;
			}
		}
		
		int score = matched - unmatched;			

		//printf("ALIGNL (%d, %d): %s\n", i, score, leftsub.c_str());
		//printf("ALIGNR (%d, %d): %s\n", i, score, rightsub.c_str());
				
		if(score > bestScore)
		{
			bestScore = score;
			bestI = i;	
		}		
		
	}

	retScore = bestScore;
	return bestI;
}

//
// Get the linkage between two specific contigs
//
ContigLinkage Scaffold::GenerateLinkage(ContigID contigID0, ContigID contigID1, PairAlignVec& paVec)
{
	ContigLinkage link;
	Contig& contig0 = m_contigMap.find(contigID0)->second;
	Contig& contig1 = m_contigMap.find(contigID1)->second;	
	assert(!paVec.empty());
	if(paVec.size() > 0)
	{
		// TODO: Check if the iter that comes back is end()
			
		// Determine the orientation of the contigs
		ContigOrientation orientation = DetermineOrientation(paVec);
		if(orientation == CORIEN_OPP)
		{
			ReverseSecondContigPairs(paVec, contig1.seq.length());
		}
		
		// Determine the order of the contigs
		ContigOrder order = DetermineOrder(paVec);
		
		// Estimate the distance between the contigs
		int distance = EstimateDistanceBetweenContigs(paVec, order, contig0.seq, contig1.seq);
		
		if(distance == -50)
		{
			link.noLink = true;
			return link;
		}
		
		link.masterID = contigID0;
		link.slaveID = contigID1;
		link.order = order;
		link.orientation = orientation;
		link.distance = distance;
		link.noLink = false;
		link.numPairs = paVec.size();
		
		if(link.numPairs > STRONG_LINK_CUTOFF)
		{
			link.type = LT_STRONG;
		}
		else
		{
			link.type = LT_WEAK;
		}
	}
	else
	{
		printf("No pairs\n");
		link.noLink = true;	
	}
	
	return link;
}

//
// Check if there are any links that are inconsistent with the chosen best link
// This would signal a repetitive, non-mergable region
//
bool Scaffold::CheckConsistency(ContigLinkage bestLink, LinkVec& alllinks)
{
	// Get the position range of the best link
	const double NUM_SIGMA = 3;
	double estDev = m_stats.GetStdDevOfEstimate(bestLink.numPairs);
	
	// Conservatively estimate the maximum position of the best link
	range bestEstimatePos;
	bestEstimatePos.start = bestLink.distance + static_cast<int>(NUM_SIGMA * estDev);
	bestEstimatePos.end = bestEstimatePos.start + m_contigMap[bestLink.slaveID].seq.length();
		
	const int CUTOFF_PAIRS = 5;
	
	// For all the links that are well supported (above the threshold cutoff), check if they substationally overlap the best link
	for(LinkIter iter = alllinks.begin(); iter != alllinks.end(); iter++)
	{
		// Ignore the self link and any poor links
		if(iter->numPairs >= CUTOFF_PAIRS && iter->slaveID != bestLink.slaveID)
		{
			// what is the minimum starting point of the current contig
			range testPosition;
			double testDev =  m_stats.GetStdDevOfEstimate(iter->numPairs);
			testPosition.start = iter->distance - static_cast<int>(NUM_SIGMA * testDev);
			testPosition.end = testPosition.start + m_contigMap[iter->slaveID].seq.length();
			
			int overlap = OverlapRanges(bestEstimatePos, testPosition);
			printf("	overlap between %s and %s is %d [%d %d] [%d %d]\n", bestLink.slaveID.c_str(), iter->slaveID.c_str(), overlap, bestEstimatePos.start, bestEstimatePos.end, testPosition.start, testPosition.end);
			if(overlap > m_kmer)
			{
				return false;
			} 
		}	
	}
	
	return true;
}

//
// Estimate the distance between the two contigs
// Precondition the contigPairs vector has the pairs in the correct (same/opp) orientation
//
int Scaffold::EstimateDistanceBetweenContigs(PairAlignVec& contigPairs, ContigOrder order, Sequence& contig1, Sequence& contig2)
{
	std::vector<int> pairDistances;
	
	int offset0 = 0;
	int offset1 = 0;
	
	if(order == CORDER_LEFT)
	{
		// Contig1 is on the left, offset contig2 by the length of contig1
		offset1 = contig1.length();
	}
	else
	{
		// Contig2 is on the left, offset contig1 by the length of contig2
		offset0 = contig2.length();
	}
	
	
	for(PAVIter iter = contigPairs.begin(); iter != contigPairs.end(); iter++)
	{
		if(iter->pairs[0].isRC != iter->pairs[1].isRC)
		{
			int distance;
			int transPos0 = iter->pairs[0].pos + offset0;
			int transPos1 = iter->pairs[1].pos + offset1;
			
			if(transPos0 < transPos1)
			{
				distance = 	transPos1 - transPos0;
			}
			else
			{
				distance = 	transPos0 - transPos1;
			}
			pairDistances.push_back(distance);
		}
	}
	
	// Perform a maximum likelihood estimate over the distances
	int estDist = m_stats.MaxLikelihoodEst(pairDistances);
	return estDist;
	
}

// 
// Assemble the vector of sequences
// 
SeqVec Scaffold::SubAssemble(PSequenceVector& seqs, Sequence startNode, Sequence stopNode, int maxDistance)
{
	// Load the phase space
	SequenceCollectionHash* pSC = new SequenceCollectionHash();
	
	// Add the reads to the phase space
	assert((int)startNode.length() == SUB_ASSEMBLY_K);
	
	
	// Add the starting node
	pSC->add(startNode);
	
	for(PSequenceVectorIterator iter = seqs.begin(); iter != seqs.end(); iter++)
	{
		// Break into kmers
		for(int i = 0; i < (m_readLen - SUB_ASSEMBLY_K + 1); i++)
		{
			PackedSeq kmer = iter->subseq(i, SUB_ASSEMBLY_K);
			pSC->add(kmer);
		}
	}	
	
	printf("Sub assembly loaded %zu pairs (max distance %d)\n", seqs.size(), maxDistance);
	printf("Stop sequence is: %s\n", stopNode.c_str());
	
	pSC->finalize();
	
	// Generate the adjacency between the sequences
	AssemblyAlgorithms::generateAdjacency(pSC);
	
	extDirection dir = SENSE;

	SeqVec assemblyProducts = AssembleRecursive(pSC, dir, PackedSeq(startNode), PackedSeq(stopNode), 0, maxDistance);
	printf("assembly returned %zu sequences\n", assemblyProducts.size());
	
	for(SeqVecIter iter = assemblyProducts.begin(); iter != assemblyProducts.end(); iter++)
	{
		printf("sub assembly returned (%zu) %s\n", iter->length(), iter->c_str());
	}
	
	delete pSC;

	return assemblyProducts;
}

//
// Recursively assemble
//
SeqVec Scaffold::AssembleRecursive(ISequenceCollection* pSC, extDirection dir, PackedSeq start, PackedSeq stop, int d, int maxDistance)
{
	PackedSeq currSeq = start;
	Sequence extSeq;
	
	while(true)
	{
		HitRecord hr = AssemblyAlgorithms::calculateExtension(pSC, currSeq, dir);
		if(hr.getNumHits() == 0 || d > maxDistance)
		{
			SeqVec ret;
			return ret;	
		}
		else if(currSeq == stop)
		{
			printf("STOP SEQ FOUND\n");
			SeqVec ret;
			ret.push_back(extSeq);
			return ret;					
		}
		else if(hr.getNumHits() == 1)
		{
			// add the extended base to the sequence
			extSeq.append(1, hr.getFirstHit().seq.getLastBase());
			currSeq = hr.getFirstHit().seq;
			++d;
		}
		else if(hr.getNumHits() > 1)
		{
			++d;
			SeqVec ret;
			
			for(int i = 0; i < hr.getNumHits(); ++i)
			{
				// TODO; Make tail recursive
				//printf("recurring at: %d\n", d);
				SeqVec outSeqs = AssembleRecursive(pSC, dir, hr.getHit(i).seq, stop, d, maxDistance);
				
				// append the current base to the extension seq
				Sequence temp = extSeq;
				temp.append(1, hr.getHit(i).seq.getLastBase());
				
				if(outSeqs.empty())
				{
					//ret.push_back(temp);
				}
				else
				{
					for(SeqVecIter iter = outSeqs.begin(); iter != outSeqs.end(); iter++)
					{
						// build the return sequence 
						ret.push_back(temp + *iter);
					}					
				}
			}
			return ret;
		}
	}
}

//
// Get all pairs from the specified contig to any other contig
//
void Scaffold::GenerateAllPairAlignments(ContigID contigID, ContigPairVecMap& cpvMap)
{
	// Get all the reads on this contig
	const ReadSet& contigReadSet = m_contigReadMap.find(contigID)->second;
	
	// Get all the pairs of this read
	for(ConstRSIter rsIter = contigReadSet.begin(); rsIter != contigReadSet.end(); rsIter++)
	{
		// Generate the pair alignment
		PairAlignVec pairAlignVec;
		bool usable = GetPairAlign(*rsIter, contigID, pairAlignVec);
		
		// Is the pair unique and both ends are aligned?
		if(usable)
		{	
			for(PAVIter pavIter = pairAlignVec.begin(); pavIter != pairAlignVec.end(); pavIter++)
			{
				if(pavIter->pairs[0].contig != pavIter->pairs[1].contig)
				{						
					// Add this pair to the vector 
					cpvMap[pavIter->pairs[1].contig].push_back(*pavIter);
				}
			}
		}
	}
	return;
	
}

//
// Get the unique pairs between the passed in contig and any other contig
//
void Scaffold::GenerateUniquePairAlignments(ContigID contigID, ContigPairVecMap& cpvMap)
{
	// Get all the reads on this contig
	const ReadSet& contigReadSet = m_contigReadMap.find(contigID)->second;
	
	// Get all the pairs of this read
	for(ConstRSIter rsIter = contigReadSet.begin(); rsIter != contigReadSet.end(); rsIter++)
	{
		// Generate the pair alignment
		PairAlign pairAlign;
		bool usable = GetUniquePairAlign(*rsIter, pairAlign);
		
		// Is the pair unique and both ends are aligned?
		if(usable)
		{	
			if(pairAlign.pairs[0].contig != pairAlign.pairs[1].contig)
			{	
				assert(pairAlign.pairs[1].contig != "");
				
				// Add this pair to the vector 
				cpvMap[pairAlign.pairs[1].contig].push_back(pairAlign);
			}
		}
	}
	return;
	
}

//
// Get all the reads of the pairs that are of the specified complement
//
void Scaffold::GetEndPairs(ContigID contigID, bool rcPairs, PSequenceVector& outSeqs)
{
	const ReadSet contigReadSet = m_contigReadMap.find(contigID)->second;

	// Get all the pairs of this read
	for(ConstRSIter rsIter = contigReadSet.begin(); rsIter != contigReadSet.end(); rsIter++)
	{
		ReadID readID = *rsIter;
		ReadID pairID = GetPairID(readID);
		
		AlignVec readAligns = GetAlignmentsForRead(readID);
		AlignVec pairAligns = GetAlignmentsForRead(pairID);
		
		assert(!readAligns.empty());
		
		// skip reads that are ambigiously mapped to this contig
		if(readAligns.size() > 1)
		{
			continue;
		}
		
		ReadAlign refAlign = readAligns.front();
		
		// Get the alignment of the pair, if its not on the same contig and its facing in the correct direction, add it to the vector
		if(!pairAligns.empty())
		{
			ReadAlign pairAlign = pairAligns.front();
			if(refAlign.contig != pairAlign.contig)
	    	{
			    // Check the pair is facing in the correct direction
			    if(refAlign.isRC == rcPairs)
			    {
					//printf("pair found on contig %s (%s)\n", pairAlign.pairs[1].contig.c_str(), pairAlign.pairs[1].seq.c_str());
					PackedSeq seq = GetSequenceForRead(pairID);
					outSeqs.push_back(seq);

				}
	    	}				
		}

    }  
}

//
//
//
bool Scaffold::GetPairAlign(ReadID readID, ContigID readContig, PairAlignVec& pairAlignVec)
{
	ReadID pairID = GetPairID(readID);
	
	AlignVec readAligns = GetAlignmentsForRead(readID);
	AlignVec pairAligns = GetAlignmentsForRead(pairID);
	
	assert(!readAligns.empty());
	
	// Ensure the root read only has 1 alignment on this contig
	AlignVec validAligns;
	for(AlignVec::iterator alignIter = readAligns.begin(); alignIter != readAligns.end(); alignIter++)
	{
		if(alignIter->contig == readContig)
		{
			validAligns.push_back(*alignIter);
		}
	}
	
	
	// Dont use ambigious pairs
	if(validAligns.empty() || validAligns.size() > 1 || pairAligns.empty())
	{
		return false;
	}
	else
	{
		ReadAlign read0Align = validAligns.front();
		
		// Create all the pairs
		for(AlignVec::iterator pairIter = pairAligns.begin(); pairIter != pairAligns.end(); pairIter++)
		{
			PairAlign temp;
			temp.pairs[0] = read0Align;
			temp.pairs[1] = *pairIter;
			temp.invalid = false;
			pairAlignVec.push_back(temp);
		}
		
		
	}	

	return true;	
}

//
// Populate the pairAlign data structure if the pairs are unique and both are aligned
// Returns true if unique/aligned, false otherwise
//
bool Scaffold::GetUniquePairAlign(ReadID readID, PairAlign& pairAlign)
{
	ReadID pairID = GetPairID(readID);
	
	AlignVec readAligns = GetAlignmentsForRead(readID);
	AlignVec pairAligns = GetAlignmentsForRead(pairID);
	
	assert(!readAligns.empty());
	
	// Dont use ambigious pairs
	if(readAligns.empty() || pairAligns.empty() || readAligns.size() > 1 || pairAligns.size() > 1)
	{
		pairAlign.invalid = true;
		return false;
	}
	else
	{
		// Create the pair align structure
		pairAlign.pairs[0] = readAligns.front();
		pairAlign.pairs[1] = pairAligns.front();
		pairAlign.invalid = false;
		return true;
	}
}

//
// Determine the orientation of the contigs
// ---------->    <------------ is CORIEN_OPP
// ---------->    ------------> is CORIEN_SAME
//
ContigOrientation Scaffold::DetermineOrientation(PairAlignVec& contigPairs)
{
	int numSameOrient = 0;
	int numOppOrient = 0;
	for(PAVIter iter = contigPairs.begin(); iter != contigPairs.end(); iter++)
	{
		if(iter->pairs[0].isRC == iter->pairs[1].isRC)
		{
			numSameOrient++;
		}
		else
		{
			numOppOrient++;
		}
	}
	
	if(numSameOrient > numOppOrient)
	{
		// The pairs are expected to be opposite orientation so if they are in the same orientation then the contigs are in the wrong orientation
		return CORIEN_OPP;
	}
	else
	{
		return CORIEN_SAME;	
	}
}

//
// Determine the order of the contigs
// The contig with more pairs in the same direction as the contig is on the LEFT since the pairs are oriented like this:
// -------->         <----------
//
ContigOrder Scaffold::DetermineOrder(PairAlignVec& contigPairs)
{
	int c1SameDir = 0;
	int c2SameDir = 0;
	
	for(PAVIter iter = contigPairs.begin(); iter != contigPairs.end(); iter++)
	{
		if(!iter->pairs[0].isRC)
		{
			c1SameDir++;
		}
		
		if(!iter->pairs[1].isRC)
		{
			c2SameDir++;
		}			
	}
	
	// The contig with more 5'->3' reads is the left contig
	if(c1SameDir > c2SameDir)
	{
		return CORDER_LEFT;
	}
	else
	{
		return CORDER_RIGHT;
	}
}

//
// Update the position of reads on the master contig
//
void Scaffold::UpdateMasterReads(ContigID contigID, int offset,
		const Sequence& /*origSeq*/, const Sequence& /*merged*/)
{
	// Get all the reads for this contig
	ReadSet& reads = m_contigReadMap[contigID];
	
	for(RSIter iter =  reads.begin(); iter != reads.end(); iter++)
	{
		//Sequence o = origSeq.substr(m_alignMap[*iter].pos, m_readLen);
		AlignVec& aligns = m_alignMap[*iter];
		for(AVIter alignIter = aligns.begin(); alignIter != aligns.end(); alignIter++)
		{
			if(alignIter->contig == contigID)
			{
				alignIter->pos += offset;
			}
		}
		
		//Sequence n = merged.substr(m_alignMap[*iter].pos, m_readLen);
		//if(o != n)
		//{
			//printf("%s == %s\n", o.c_str(), n.c_str());
			//assert(o == n);
		//}
	}
}

//
// Update the reads of the contig that got merged into the master
//
void Scaffold::UpdateSlaveReads(ContigID slaveID, ContigID masterID,
		int offset, bool isFlipped, const Sequence& origSeq,
		const Sequence& /*merged*/)
{
	// Get all the reads for this contig
	ReadSet& reads = m_contigReadMap[slaveID];
	
	for(RSIter iter =  reads.begin(); iter != reads.end(); iter++)
	{
		//Sequence o = origSeq.substr(m_alignMap[*iter].pos, m_readLen);
		AlignVec& aligns = m_alignMap[*iter];
		for(AVIter alignIter = aligns.begin(); alignIter != aligns.end(); alignIter++)
		{		
			int origPos = alignIter->pos;
			if(isFlipped)
			{
				// flip the reads position 
				alignIter->pos = origSeq.length() - origPos - m_readLen;
				alignIter->isRC = !alignIter->isRC;
			}
			
			// Offset the position
			alignIter->pos += offset;
			alignIter->contig = masterID;
			
			// Add the read id to the master contig
			m_contigReadMap[masterID].insert(*iter);
		}
		
		//Sequence n = merged.substr(m_alignMap[*iter].pos, m_readLen);
		//printf("%s == %s (%d %d %d)\n", o.c_str(), n.c_str(), offset, isFlipped, origPos);
		//assert(o == n);
	}
	
	// Remove the slave contig's alignments
	m_contigReadMap.erase(slaveID);
	
}

//
// Realign the pairs of the reads on the specified contig
//
void Scaffold::RealignContigPairs(ContigID /*contigID*/)
{
	assert(false);
	/*
	// Create the aligner
	printf("Creating DB\n");	
	Aligner aligner(29);
	ContigMap tempMap;
	tempMap[contigID].seq = m_contigMap[contigID].seq;
	aligner.CreateDatabase(tempMap);	
	
	// Get all the reads for this contig
	ReadSet& reads = m_contigReadMap[contigID];
	
	// Iterate over all the reads in the set
	for(RSIter iter =  reads.begin(); iter != reads.end(); iter++)
	{	
		// Get the reads and their alignments
		ReadID readID = *iter;
		ReadID pairID = GetPairID(readID);
		
		AlignVec readAligns = GetAlignmentsForRead(readID);
		AlignVec pairAligns = GetAlignmentsForRead(pairID);
		
		// If the pair is uniquely aligned to this contig, nothing has to be done
		if(pairAligns.size() == 1 && pairAligns.front().contig == contigID)
		{
			continue;	
		}
		else
		{		
			PackedSeq seq = GetSequenceForRead(pairID);
			
			AlignmentVector av = aligner.GetAlignments(PackedSeq(seq));
			
			if(!av.empty())
			{
				// This pair is mapped to the new contig, delete its old alignments
				for(AVIter oldAlignIter = pairAligns.begin(); oldAlignIter != pairAligns.end(); oldAlignIter++)
				{
					m_contigReadMap[oldAlignIter->contig].erase(pairID);
				}
				
				// Delete the old alignments record
				m_alignMap.erase(pairID);
				
				// Create the new alignment(s)
				for(AlignmentVector::iterator arIter = av.begin(); arIter != av.end(); arIter++)
				{
					ReadAlign ra = BuildReadAlign(pairID, arIter->contig, arIter->start, arIter->isRC);
					
					// Set up all the mappings
					m_alignMap[pairID].push_back(ra);
					
					m_contigReadMap[arIter->contig].insert(pairID);
				}
			}
		}
	}
	*/
}

//
// Reverse the pairs of the second contig to make the orientations the same
//
void Scaffold::ReverseSecondContigPairs(PairAlignVec& contigPairs, int contigLength)
{
	for(PAVIter iter = contigPairs.begin(); iter != contigPairs.end(); iter++)
	{
		iter->pairs[1].pos = contigLength - iter->pairs[1].pos - m_readLen;
		iter->pairs[1].isRC = !iter->pairs[1].isRC;
	}
}	

//
//
//
int Scaffold::OverlapRanges(const range r1, const range r2)
{
	if(r1.start < r2.start)
	{
		return r1.end - r2.start;
	}
	else
	{
		return r2.end - r1.start;
	}
}

//
//
//
range Scaffold::GenerateRange(int distance, int size, int n, int numDevs)
{
	double s = m_stats.GetStdDevOfEstimate(n);
	int error = (int)s * numDevs;
	range r;
	r.start = distance - error;
	r.end = distance + error + size;
	return r;
}

//
//
//
AlignVec Scaffold::GetAlignmentsForRead(ReadID id)
{
	PairAlign pairAlign;
	
	// Get the read alignments
	ConstAMIter raIter = m_alignMap.find(id);
	if(raIter != m_alignMap.end())
	{
		return raIter->second;
	}
	else
	{
		AlignVec empty;
		return empty;	
	}
}

//
// Get the sequence for a particular read
//
PackedSeq Scaffold::GetSequenceForRead(ReadID id)
{
	// Since the ID is simply an int and the reads are loaded in order of their id, simply return the sequence at index id
	assert(id >= 0 && id < (int)m_readVec.size());
	return m_readVec[id];
}

//
// Get the ID of a pair from a readID
//
ReadID Scaffold::GetPairID(ReadID id)
{
	if(id % 2 == 0)
	{
		// even id
		return (id + 1);
	}
	else
	{
		// odd id
		return (id - 1);
	}
}

//
// Build the empirical distribution of pair distances from contigs that have both ends of the pair in the contig
//
void Scaffold::GenerateStatistics()
{
	
	std::map<int, int> distribution;
	for(ConstCRMIter iter = m_contigReadMap.begin(); iter != m_contigReadMap.end(); iter++)
	{
		const ContigID id = iter->first;
		
		// Get the vector of reads for this contig
		const ReadSet& readSet = iter->second;
		
		// Iterate over the reads, checking which ones have a pair on the same contig
		for(ConstRSIter iter2 = readSet.begin(); iter2 != readSet.end(); iter2++)
		{
			// Generate the pair alignment
			PairAlign pairAlign;
			bool usable = GetUniquePairAlign(*iter2, pairAlign);
			if(usable)
			{
				if(pairAlign.pairs[0].contig == pairAlign.pairs[1].contig && pairAlign.pairs[0].isRC != pairAlign.pairs[1].isRC)
				{
					int distance;
					if(pairAlign.pairs[0].pos < pairAlign.pairs[1].pos)
					{
						distance = pairAlign.pairs[1].pos - pairAlign.pairs[0].pos;
					}
					else
					{
						distance = pairAlign.pairs[0].pos - pairAlign.pairs[1].pos;
					}
					distribution[distance]++;
				}
			}
		}
	}
	
	// set up the stats
	m_stats.GenerateStatsFromHistogram(distribution);
}

//
// Write the alignments structure to a file
//
void Scaffold::WriteAlignments(std::string filename)
{
	std::ofstream fileHandle(filename.c_str());	
	for(AMIter readIter = m_alignMap.begin(); readIter != m_alignMap.end(); readIter++)
	{
		for(AlignVec::iterator alignIter = readIter->second.begin(); alignIter != readIter->second.end(); alignIter++)
		{
			fileHandle << readIter->first << " " << alignIter->contig << " " << alignIter->pos << " " << alignIter->isRC << " " << "NOSEQ" << endl;
		}
	}
	fileHandle.close();
	//m_alignMap[id].push_back(ra);
}

// Load the adjacency info
void Scaffold::LoadAdjacency(std::string file)
{
	m_pSC = new SequenceCollectionHash();
	
	AssemblyAlgorithms::loadSequences(m_pSC, file);
	m_pSC->finalize();
	AssemblyAlgorithms::generateAdjacency(m_pSC);
	
	printf("Done reconstructing DBG\n");
}

//
// Read in a fasta file of sequence reads
//
void Scaffold::ReadSequenceReads(std::string file)
{
	// Read in the fasta file
	FastaReader reader(file.c_str());
	int count = 0;
	while(reader.isGood())
	{
		PackedSeq seq = reader.ReadSequence();
		m_readVec.push_back(seq);
		count++;
	}
	printf("Read %d sequences\n", count);
}

//
// Read a pairing file
//
void Scaffold::ReadPairs(std::string file)
{
	// Read in the pairs file
	std::ifstream fileHandle(file.c_str());	
	while(!(fileHandle.eof() || fileHandle.peek() == EOF))
	{
		//fileHandle.getline(buffer, MAX_LINE_LEN);	
		
		ReadID id1;
		ReadID id2;
		std::string c1;
		std::string c2;
		int pos1;
		int pos2;
		int rc1;
		int rc2;
		Sequence seq1;
		Sequence seq2;
		
		fileHandle >> id1 >> c1 >> pos1 >> rc1 >> seq1 >> id2 >> c2 >> pos2 >> rc2 >> seq2;

		// Check if this read was valid
		if(fileHandle.eof())
		{
			continue;
		}
		
		// don't record the sequences for pairs on the same chromosomes (a bit hacky of a way to save memory)
		if(c1 == c2)
		{
			seq1 = "";
			seq2 = "";
		}
		
		ReadAlign ra1 = BuildReadAlign(id1, c1, pos1, rc1);
		ReadAlign ra2 = BuildReadAlign(id2, c2, pos2, rc2);
		
		// Set up all the mappings
		m_alignMap[id1].push_back(ra1);
		m_alignMap[id2].push_back(ra2);
		
		m_contigReadMap[c1].insert(id1);
		m_contigReadMap[c2].insert(id2);
	}
}


void Scaffold::ReadAlignments(std::string file)
{
	// Read in the alignments file
	std::ifstream fileHandle(file.c_str(), ios_base::in);
	
	assert(fileHandle.is_open());
	
	while(fileHandle.good())
	{
		//fileHandle.getline(buffer, MAX_LINE_LEN);	
		ReadID id;
		std::string c;
		int pos;
		int rc;
		Sequence seq;

		fileHandle >> id >> c >> pos >> rc >> seq;

		// Check if this read was valid
		if(fileHandle.eof())
		{
			continue;
		}

		ReadAlign ra = BuildReadAlign(id, c, pos, rc);
		
		// Set up all the mappings
		m_alignMap[id].push_back(ra);
		
		m_contigReadMap[c].insert(id);
	}
	
	fileHandle.close();
}

//
// Build a readalign structure
//

ReadAlign Scaffold::BuildReadAlign(ReadID id, std::string contig, int position, bool isRC)
{
	ReadAlign ra;
	ra.id = id;
	ra.contig = contig;
	ra.pos = position;
	ra.isRC = isRC;	
	return ra;
}

//
// Print a read alignment
//
void Scaffold::PrintReadAlign(ReadAlign& ra)
{
	std::cout << "Read: " << ra.id << " c: " << ra.contig << " p: " << ra.pos << " rc: " << ra.isRC << std::endl;
}

// 
int CompareLinkagesByDistance(const ContigLinkage& l1, const ContigLinkage& l2)
{
	return l1.distance < l2.distance;
}

// 
int CompareLinkagesByDistanceDesc(const ContigLinkage& l1, const ContigLinkage& l2)
{
	return l1.distance > l2.distance;
}


//
// Generate the graphviz-readable graph file for the contig
//
void Scaffold::GenerateGraph(ContigID contigID)
{
	ContigPairVecMap cpvMap;
	// Get all the pair alignments for this contig
	// It will generate a mapping of contig pairs (c1,c2) to a vector of the pairs supporting the contigs
	// This function populates the cpvMap structure
	GenerateUniquePairAlignments(contigID, cpvMap);
	
	ofstream file("contigGraph.dot");
	file << "digraph G { size = \"8.5, 11.0!\"; ratio = \"fill\";" << endl;
	
	std::set<std::string> seenLinks;
	
	for(CPVMIter cpvmIter = cpvMap.begin(); cpvmIter != cpvMap.end(); cpvmIter++)
	{		
		ContigLinkage link = GenerateLinkage(contigID, cpvmIter->first, cpvmIter->second);
		
		std::string orderString;
		
		if(link.order != CORDER_RIGHT)
		{
			continue;
		}
		OutputGVizNode(file, link);
		//continue;
		ContigPairVecMap cpvMap2;
		GenerateUniquePairAlignments(link.slaveID, cpvMap2);
		
		for(CPVMIter cpvmIter2 = cpvMap2.begin(); cpvmIter2 != cpvMap2.end(); cpvmIter2++)
		{
			if(cpvmIter2->first != contigID)
			{
				std::string seenstr;
				if(link.slaveID < cpvmIter2->first)
				{
					seenstr = link.slaveID;
					seenstr += "--";
					seenstr += cpvmIter2->first;
				}
				
				// avoid duplicates
				if(seenLinks.find(seenstr) == seenLinks.end())
				{
					ContigLinkage link2 = GenerateLinkage(link.slaveID, cpvmIter2->first, cpvmIter2->second);
					OutputGVizNode(file, link2);
					seenLinks.insert(seenstr);
				}
			}
		}
	}

	
	file << "}" << endl;
}

void Scaffold::OutputGVizNode(std::ofstream& ostr, ContigLinkage& link)
{
	std::string orderString;
	
	if(link.order == CORDER_LEFT)
	{
		orderString = link.masterID +  " -> " +  link.slaveID;
	}
	else
	{
		orderString = link.slaveID +  " -> " + link.masterID;
	}	
	
	range r = GenerateRange(link.distance, m_contigMap[link.slaveID].seq.length(), link.numPairs, 1);
	ostr << "\t" << orderString << " [arrowhead = \"normal\" label = \"[" << r.start << "," << r.end << "], " << link.numPairs << "\"];" << endl;
}

void Scaffold::LoadPairsRecord(const PSequenceVector& allreads, int kmerSize)
{
	size_t numReads = allreads.size();
	int count = 0;
	for(size_t index = 0; index < numReads; index += 2)
	{
		Sequence seq1 = allreads[index].decode();
		Sequence seq2 = allreads[index+1].decode();
				
		int len = seq1.length();
		assert(kmerSize <= len);
		
		for(int i = 0; i < len - kmerSize  + 1; i++)
		{
			PackedSeq sub1 = seq1.substr(i, kmerSize);
			PackedSeq sub2 = seq2.substr(i, kmerSize);
			
			m_pairRec.addPairs(sub1, sub2);
			count++;
		}
	}
	printf("Read %d pairs\n", count);
}

//
// Print a linkage object
//
void Scaffold::PrintLinkage(ContigLinkage& link)
{
	std::string strOrient = (link.orientation == CORIEN_SAME) ? "same" : "opposite";
	std::string strOrder = (link.order == CORDER_LEFT) ? "left" : "right";
	cout << '\t' << link.slaveID << " orientation: " << strOrient << " c1 is on " << strOrder << " distance: " << link.distance << endl;	
}
