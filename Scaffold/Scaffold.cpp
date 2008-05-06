#include <stdio.h>
#include <math.h>
#include "Scaffold.h"
#include "Sequence.h"
#include "FastaWriter.h"
#include "SequenceCollectionHash.h"
#include "Aligner.h"
#include "FastaReader.h"

using namespace std;

int main(int argc, char** argv)
{
	if(argc < 4)
	{
		printf("usage: Scaffold <readsFile> <contigFile> <readLen> <assembly kmer size> [alignment file]\n");
		exit(1);
	}
	
	std::string readsFile(argv[1]);
	std::string contigFile(argv[2]);
	int readLen = atoi(argv[3]);
	int kmer = atoi(argv[4]);
	
	// Should we read in the alignment hints file?
	bool readAlignments = false;
	std::string alignmentsFile = "";
	if(argc > 5)
	{
		readAlignments = true;
		alignmentsFile = argv[5];
	}
	
	Scaffold scaffold(readsFile, contigFile, readLen, kmer, readAlignments, alignmentsFile);
}

//
//
//
Scaffold::Scaffold(std::string readsFile, std::string contigFile, int readLen, int kmer, bool bReadAlignments, std::string alignmentsFile) : m_readLen(readLen), m_kmer(kmer)
{
	// Read in the pairs
	//ReadAlignments(alignFile);
	
	// Read in the sequencing reads
	ReadSequenceReads(readsFile);
	
	// Read in the contigs
	ReadContigs(contigFile);	
	
	// Generate the initial alignment database
	
	if(!bReadAlignments)
	{
		printf("Generating alignments from reads file\n");
		
		GenerateAlignments(m_readVec, m_contigMap);
		WriteAlignments("alignments.txt");
	}
	else
	{
		printf("Reading the alignments from file %s\n", alignmentsFile.c_str());
		ReadAlignments(alignmentsFile);
	}
	

	
	printf("Done reading\n");
	// Generate the empirical distribution
	
	GenerateStatistics();

	while(AttemptMerge("731"));
	exit(1);
	//AttemptMerge("1219");
	//exit(1);
	/*
	AttemptMerge("1128");
	AttemptMerge("2523");
	AttemptMerge("1476");
	AttemptMerge("2823");
	//GenerateGraph("1009");	
	exit(1);
	*/
	
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

//
// Generate the initial alignments for the sequences
//
void Scaffold::GenerateAlignments(PSequenceVector& seqs, ContigMap& /*contigs*/)
{
	printf("Creating DB\n");	
	Aligner aligner(29);
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
	
				// Get all the pairs for these contigs that fall into the gap between the contigs
				PSequenceVector pairs;
				GetEndPairs(bestLink->masterID, contig0Comp, pairs);
				GetEndPairs(bestLink->slaveID, contig1Comp, pairs);
				
				// Now the contigs are in the correct order (leftContig->rightContig) and are from the same strand
				
				// Attempt to assemble the gap between the contigs using the gathered pair info
					
				// generate the sequence to start from
				// this is the last K bases of the reference (contig0) sequence
				Sequence start = leftContig->substr(leftContig->length() - SUB_ASSEMBLY_K, SUB_ASSEMBLY_K);
				
				// generate the stop sequence (the first k bases of the right contig)
				Sequence stop = rightContig->substr(0, SUB_ASSEMBLY_K);
				
				printf("Starting assembly from %s\n", start.c_str());
				
				// Assemble the sequence
				int expectedAssemblyLength = bestLink->distance + SUB_ASSEMBLY_K;
				int range = 20;
				SeqVec assemblies = SubAssemble(pairs, start, stop, expectedAssemblyLength + range);
					
				// Check if any of the returned assemblies are valid
				int bestScore = abs(bestLink->distance);
				SeqVecIter bestIter;
				bool hasAssembly = false;
				
				// Choose the assembly that is closest in size to the expected distance
				for(SeqVecIter iter = assemblies.begin(); iter != assemblies.end(); iter++)
				{	
					// Estimate the position of the assembly on the right contig
					int len = iter->length();
					int score = abs(len - expectedAssemblyLength);
					printf("score: %d best: %d\n", score, bestScore);
					if(score < bestScore)
					{
						bestScore = score;
						bestIter = iter;
						hasAssembly = true;
					}
				}
				
				Sequence merged;
				if(hasAssembly)
				{
					// Merge with the generated path
					Sequence tempMerged;
					
					// Merge the left contig and the assembly first at distance = 0
					Merge(*leftContig, *bestIter, 0, tempMerged);
					
					// The assembly product will overlap the right contig by SUB_ASSEMBLY_K by definition
					int mergeDistance = -SUB_ASSEMBLY_K;
					
					// Merge the right contig
					Merge(tempMerged, *rightContig, mergeDistance, merged);
				}
				else
				{
					// Merge with a gap	
					Merge(*leftContig, *rightContig, bestLink->distance, merged);
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
			}
			
			printf("MERGED (%zu): %s\n", m_contigMap[contigID].seq.length(), m_contigMap[contigID].seq.c_str());
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
void Scaffold::RealignContigPairs(ContigID contigID)
{
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
// Read contigs
//
void Scaffold::ReadContigs(std::string file)
{
	std::ifstream fileHandle(file.c_str());	
	while(fileHandle.good())
	{
		char head;

		std::string contigID;
		//int size;
		//double coverage;
		Sequence seq;
		fileHandle >> head;
		
		if(head != '>' || fileHandle.eof())
		{
			continue;
		}		
		
		
		fileHandle >> contigID;
		fileHandle.ignore(1000, '\n');	
		fileHandle >> seq;

		m_contigMap[contigID].seq = seq;
		m_contigMap[contigID].merged = false;
		m_contigMap[contigID].repetitive = false;
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

//
// Print a linkage object
//
void Scaffold::PrintLinkage(ContigLinkage& link)
{
	std::string strOrient = (link.orientation == CORIEN_SAME) ? "same" : "opposite";
	std::string strOrder = (link.order == CORDER_LEFT) ? "left" : "right";
	cout << '\t' << link.slaveID << " orientation: " << strOrient << " c1 is on " << strOrder << " distance: " << link.distance << endl;	
}
