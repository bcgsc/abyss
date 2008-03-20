#include <stdio.h>
#include <math.h>
#include "Scaffold.h"
#include "Sequence.h"
#include "FastaWriter.h"

using namespace std;

const double MINP = 0.00001f;

int main(int argc, char** argv)
{
	std::string pairsFile(argv[1]);
	std::string contigFile(argv[2]);
	int readLen = atoi(argv[3]);
	int kmer = atoi(argv[4]);
	
	Scaffold scaffold(pairsFile, contigFile, readLen, kmer);
}

//
//
//
Scaffold::Scaffold(std::string pairsFile, std::string contigFile, int readLen, int kmer) : m_readLen(readLen), m_kmer(kmer)
{
	// Read in the pairs
	ReadPairs(pairsFile);
	
	// Read in the contigs
	ReadContigs(contigFile);	
	
	printf("Done reading\n");
	// Generate the empirical distribution
	
	GenerateEmpDistribution();
	GenerateGraph("1855");
	exit(1);
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
			if(!iter->second.merged && iter->second.seq.length() > 300)
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
// Get all the linkages between this contig and its pairs
//
bool Scaffold::AttemptMerge(ContigID contigID)
{
	bool mergedOccured = false;
	ContigPairVecMap cpvMap;
	// Get all the pair alignments for this contig
	// It will generate a mapping of contig pairs (c1,c2) to a vector of the pairs supporting the contigs
	// This function populates the cpvMap structure
	GenerateAllPairAlignments(contigID, cpvMap);
	
	cout << "Contig " << contigID << " is linked to: " << endl;
	
	// Find the linkages of this contig on the left and right
	LinkVec linkages[2];
	
	for(CPVMIter cpvmIter = cpvMap.begin(); cpvmIter != cpvMap.end(); cpvmIter++)
	{		
		ContigLinkage link = GenerateLinkage(contigID, cpvmIter->first, cpvmIter->second);
		if(!link.noLink)
		{
			linkages[link.order].push_back(link);
		}
	}
	
	for(int i = 0; i <= 1; i++)
	{
		printf("Contig %s (%d bp) has %d linkages on the %d side\n", contigID.c_str(), m_contigMap[contigID].seq.length(), linkages[i].size(), i);

		//while(RefineLinkages(linkages[i]));
		bool consistent = true;
		int prevEnd = 0;
		
		// Sort the vector
		sort(linkages[i].begin(), linkages[i].end(), CompareLinkagesByDistance);

		for(LinkIter iter = linkages[i].begin(); iter != linkages[i].end(); iter++)
		{
			
			if(iter->noLink)
			{
				consistent = false;
				continue;
			}
			
			Sequence& contig2 = m_contigMap.find(iter->slaveID)->second.seq;			
			int start = iter->distance;
			int end = start + contig2.length();
			int overlap = start - prevEnd;
							
			printf("Contig %s is %d bp away starts at %d and ends at %d (overlap %d) (pairs %d)\n", iter->slaveID.c_str(), iter->distance, start, end, overlap, iter->numPairs);
			if(iter->type == LT_STRONG)
			{		
				if(-overlap >= (m_kmer + 10))
				{
					consistent = false;
				}
				
				prevEnd = end;
			}	
		}

		if(consistent)
		{
			// The path is consistent, merge the contigs along it
			int offsetDistance = 0;
			for(std::vector<ContigLinkage>::iterator iter = linkages[i].begin(); iter != linkages[i].end(); iter++)
			{
				
				Sequence merged;
				// Offset the distance between the contigs by the total amount that has been merged in
				iter->distance -= offsetDistance;
				
				
				if(iter->type == LT_STRONG)
				//if(iter->distance < 0 || iter->type == LT_STRONG)
				{
					printf("Merging %s in at distance %d (offset: %d)\n", iter->slaveID.c_str(), iter->distance, offsetDistance);					
					// Perform the merge and update the distance
					offsetDistance += Merge(*iter, merged);
										
					int c0Len = m_contigMap[contigID].seq.length();
					int c1Len = m_contigMap[iter->slaveID].seq.length();
					
					// Update the pairs record
					// If this was the left contig, no update needs to happen, else update by the amount of sequence gained
					if(iter->order == CORDER_RIGHT)
					{
						// Offset the master contig pairs by the increase in distance
						int offset = merged.length() - c0Len;
						UpdateMasterReads(contigID, offset, m_contigMap[contigID].seq, merged);
					}
					
					// Update the slave record
					int offset = 0;
					if(iter->order == CORDER_LEFT)
					{
						offset = merged.length() - c1Len;
					}
					
					bool isFlipped = iter->orientation == CORIEN_OPP;
					UpdateSlaveReads(iter->slaveID, contigID, offset, isFlipped, m_contigMap[iter->slaveID].seq, merged);
					
					// Replace the contig with the merged product
					m_contigMap[contigID].seq = merged;
					
					// Erase the existance of the slave contig
					m_contigMap.erase(iter->slaveID);
					
					mergedOccured = true;
				}
				
			}
		}
		
		printf("MERGED (%d): %s\n", m_contigMap[contigID].seq.length(), m_contigMap[contigID].seq.c_str());
	}
	
	return mergedOccured;
}

//
// Merge two contigs using the information provided in the link data structure
//
int Scaffold::Merge(ContigLinkage& link, Sequence& merged)
{
	Sequence& contig0 = m_contigMap.find(link.masterID)->second.seq;
	Sequence& contig1 = m_contigMap.find(link.slaveID)->second.seq;
	
	printf("length (%d %d)\n", contig0.length(), contig1.length());
	
	if(link.orientation == CORIEN_OPP)
	{
		// Swap the orientation of the slave contig
		contig1 = reverseComplement(contig1);
	}
	
	Sequence leftContig;
	Sequence rightContig;
	if(link.order == CORDER_LEFT)
	{
		// The master contig is first
		leftContig = contig0;
		rightContig = contig1;
	}
	else
	{
		leftContig = contig1;
		rightContig = contig0;	
	}
	

	// Check if the contigs (potentially overlap)

	if(link.distance < 0)
	{
		int startAlign = link.distance - 4*((int)(m_stdDev / sqrt(link.numPairs)));
		if(startAlign < -(m_kmer - 1))
		{
			startAlign = -(m_kmer - 1);
		}
		int score;
		int alignDistance = partialAlign(leftContig, rightContig, startAlign, score);
		
		// Don't trust mediocre hits
		if(score < 10)
		{
			alignDistance = link.distance;
		}

		// Generate the merged sequence
		merged = leftContig;
		
		// Add the substring of the second sequence
		merged.append(rightContig.substr(abs(alignDistance), rightContig.length() - abs(alignDistance)));

		//printf("ALIGN DIST: %d\n", alignDistance);
	}
	else
	{
		// Positive distance, fill the gap with "N"
		merged = leftContig;
		merged.append(link.distance, 'N');
		merged.append(rightContig);
	}
	
	
	int distDiff = merged.length() - contig0.length();
	
	printf("MERGED AT: %d\n", merged.length() - (contig0.length() + contig1.length()));
	return distDiff;
}

//
// Align the ends of two contigs based on the passed in distance and return a refined distance for the maximum overlap
//
int Scaffold::partialAlign(const Sequence& leftContig, const Sequence& rightContig, int startPos, int &retScore)
{
	int bestScore = -1000;
	int start = startPos;
	int bestI = start;
	//printf("Aligning: %s with %s estimated to be %d\n", leftContig.c_str(), rightContig.c_str(), estDistance);
	
	for(int i = start; i < 0; i++)
	{
		int l = abs(i);
		
		// Cut off the last l bases of leftContig and the first l bases of rightContig
		Sequence leftsub = leftContig.substr(leftContig.length() - l, l);
		Sequence rightsub = rightContig.substr(0, l);
		assert(leftsub.length() == rightsub.length());
		
		int matched = 0;
		int unmatched = 0;
		for(int j = 0; j < l; j++)
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
		
		printf("ALIGNL (%d, %d): %s\n", i, score, leftsub.c_str());
		printf("ALIGNR (%d, %d): %s\n", i, score, rightsub.c_str());
		
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
	double estDev = GetStdDevOfEstimate(bestLink.numPairs);
	
	// Conservatively estimate the maximum position of the best link
	range bestEstimatePos;
	bestEstimatePos.start = bestLink.distance + static_cast<int>(NUM_SIGMA * estDev);
	bestEstimatePos.end = bestEstimatePos.start + m_contigMap[bestLink.slaveID].seq.length();
		
	const int CUTOFF_PAIRS = 10;
	
	// For all the links that are well supported (above the threshold cutoff), check if they substationally overlap the best link
	for(LinkIter iter = alllinks.begin(); iter != alllinks.end(); iter++)
	{
		if(iter->numPairs >= CUTOFF_PAIRS)
		{
			// what is the minimum starting point of the current contig
			range testPosition;
			double testDev =  GetStdDevOfEstimate(iter->numPairs);
			testPosition.start = iter->distance - static_cast<int>(NUM_SIGMA * testDev);
			testPosition.end = testPosition.start + m_contigMap[iter->slaveID].seq.length();
			
			int overlap = OverlapRanges(bestEstimatePos, testPosition);
			if(-overlap > m_kmer)
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
	int estDist = MaxLikelihoodEst(pairDistances, m_pdf);
	return estDist;
	
}

//
// Get all the pairs between the passed in contig and any other contig
//
void Scaffold::GenerateAllPairAlignments(ContigID contigID, ContigPairVecMap& cpvMap)
{
	// Get all the reads on this contig
	const ReadVec contigReadVec = m_contigReadMap.find(contigID)->second;
	
	// Get all the pairs of this read
	for(ConstRVIter rvIter = contigReadVec.begin(); rvIter != contigReadVec.end(); rvIter++)
	{
		PairAlign pairAlign = GetPairAlignsForRead(*rvIter);
		if(pairAlign.pairs[0].contig != pairAlign.pairs[1].contig)
		{
			// Add this pair to the vector 
			cpvMap[pairAlign.pairs[1].contig].push_back(pairAlign);
		}
	}
	
	return;
	
}

//
// Get the pairs between the contigs
//
void Scaffold::GeneratePairAlignmentsBetweenContigs(ContigID contigID0, ContigID contigID1, PairAlignVec& paVec)
{
	// Get all the reads on this contig
	const ReadVec contigReadVec = m_contigReadMap.find(contigID0)->second;
	
	// Get all the pairs of this read
	for(ConstRVIter rvIter = contigReadVec.begin(); rvIter != contigReadVec.end(); rvIter++)
	{
		PairAlign pairAlign = GetPairAlignsForRead(*rvIter);
		if(pairAlign.pairs[0].contig == contigID0 && pairAlign.pairs[1].contig == contigID1)
		{
			// Add this pair to the vector 
			paVec.push_back(pairAlign);
		}
	}
	
	return;	
	
}

//
// Refine links
//
bool Scaffold::RefineLinkages(LinkVec& links)
{
	assert(false);
	bool seqRefined = false;
	for(LinkIter strIter = links.begin(); strIter != links.end(); strIter++)
	{
		if(strIter->type == LT_STRONG)
		{
			Sequence strongContig = m_contigMap[strIter->slaveID].seq;
				
			int strongStart = strIter->distance;
			int strongEnd = strongStart + strongContig.length();
			
			// Iterate over the strong links, attempting to refine the position of weak links
			for(LinkIter subIter = links.begin(); subIter != links.end(); subIter++)
			{
				//Skip if its the same link or another strong link
				if(strIter == subIter || subIter->type == LT_STRONG)
				{
					continue;
				}
				

				Sequence weakContig =  m_contigMap[subIter->slaveID].seq;				
				
				int subStart = subIter->distance;
				//int subEnd = subStart + weakContig.length();
				
				// Check if these contigs overlap within a reasonable tolerance
				//bool posOverlap = (strongStart > subStart) ? strongStart < subEnd : subStart < strongEnd; 

				if(strIter->orientation != subIter->orientation)
				{
					weakContig = reverseComplement(weakContig);
				}
				
				Sequence lcontig;
				Sequence rcontig;
				if(subStart < strongStart)
				{
					lcontig = weakContig;
					rcontig = strongContig;
				}
				else
				{
					lcontig = strongContig;
					rcontig = weakContig;	
					
				}
				

				int	startAlign = -(m_kmer - 1);			
				int score;
				int overlap = partialAlign(lcontig, rcontig, startAlign, score);
				printf("Refine: alignment between (%s, %s) produces overlap of %d\n", subIter->slaveID.c_str(), strIter->slaveID.c_str(), overlap);
				// Is the overlap sufficient?
				if(abs(overlap) > 15)
				{
					if(subStart < strongStart)
					{
						subIter->distance = strongStart - (weakContig.length() + overlap);
					}
					else
					{
						subIter->distance = strongEnd + overlap;
					}
					subIter->type = LT_STRONG;
					printf("REFINED ALIGNMENT OF %s using %s to distance: %d\n", subIter->slaveID.c_str(), strIter->slaveID.c_str(), subIter->distance);
					seqRefined = true;
				}		
			}
		}
	}
	return seqRefined;
}

int Scaffold::MaxLikelihoodEst(std::vector<int>& pairDistance, PDF& pdf)
{
	double maxL = -999999;
	int bestDist = -50;
	
	// TODO: Unhardcode this
	int minDist = -50;
	int maxDist = 200;
	for(int i = minDist; i < maxDist; i++)
	{
		double v = 	ComputeLikelihood(i, pairDistance, pdf);
		if(v > maxL)
		{
			maxL = v;
			bestDist = i;
		}
	}
	
	return bestDist;
}

//
// Compute the log likelihood function over the test distribution
//
double Scaffold::ComputeLikelihood(int d, std::vector<int>& testDist, PDF& pdf)
{
	double sum = 0.0f;

	for(vector<int>::iterator iter = testDist.begin(); iter != testDist.end(); iter++)
	{
		int val = *iter + d;
		double p;
		
		// Is this value in range of the pdf?
		if(val >= 0 && val < pdf.size())
		{
			p = pdf[val];	
		}
		else
		{
			p = MINP;
		}
		
		sum += log(p);
	}
	
	return sum;
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

void Scaffold::UpdateMasterReads(ContigID contigID, int offset, const Sequence& origSeq, const Sequence& merged)
{
	// Get all the reads for this contig
	ReadVec& reads = m_contigReadMap[contigID];
	
	for(RVIter iter =  reads.begin(); iter != reads.end(); iter++)
	{
		//Sequence o = origSeq.substr(m_alignMap[*iter].pos, m_readLen);
		m_alignMap[*iter].pos += offset;
		//Sequence n = merged.substr(m_alignMap[*iter].pos, m_readLen);
		//if(o != n)
		//{
			//printf("%s == %s\n", o.c_str(), n.c_str());
			//assert(o == n);
		//}
	}
}


void Scaffold::UpdateSlaveReads(ContigID slaveID, ContigID masterID, int offset, bool isFlipped, const Sequence& origSeq, const Sequence& merged)
{
	// Get all the reads for this contig
	ReadVec& reads = m_contigReadMap[slaveID];
	
	for(RVIter iter =  reads.begin(); iter != reads.end(); iter++)
	{
		//Sequence o = origSeq.substr(m_alignMap[*iter].pos, m_readLen);
		int origPos = m_alignMap[*iter].pos;
		if(isFlipped)
		{
			// flip the reads position 
			m_alignMap[*iter].pos = origSeq.length() - origPos - m_readLen;
			m_alignMap[*iter].isRC = !m_alignMap[*iter].isRC;
		}
		
		// Offset the position
		m_alignMap[*iter].pos += offset;
		m_alignMap[*iter].contig = masterID;
		
		// Add the read to the master contig
		m_contigReadMap[masterID].push_back(*iter);
		
		//Sequence n = merged.substr(m_alignMap[*iter].pos, m_readLen);
		//printf("%s == %s (%d %d %d)\n", o.c_str(), n.c_str(), offset, isFlipped, origPos);
		//assert(o == n);
	}
	
	// Remove the slave contig's alignments
	m_contigReadMap.erase(slaveID);
	
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
// Build the empirical distribution of pair distances from contigs that have both ends of the pair in the contig
//
void Scaffold::GenerateEmpDistribution()
{
	
	std::map<int, int> distribution;
	for(ConstCRMIter iter = m_contigReadMap.begin(); iter != m_contigReadMap.end(); iter++)
	{
		const ContigID id = iter->first;
		
		// Get the vector of reads for this contig
		const ReadVec& readVec = iter->second;
		
		// Iterate over the reads, checking which ones have a pair on the same contig
		for(ConstRVIter iter2 = readVec.begin(); iter2 != readVec.end(); iter2++)
		{
			ReadID r1ID = *iter2;
			PairAlign pairAlign = GetPairAlignsForRead(r1ID);
			
			// Ensure these reads are on the same contig and have opposite orientations
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
	
	// set m_pdf
	ConvertHistToPDF(distribution);
}

//
// Convert the histogram to a probability density function
// DOUBLE COUNTING PAIRS??
void Scaffold::ConvertHistToPDF(const histogram& h)
{
	int max = 0;
	double count = 0;
	double sum = 0;
	for(histogram::const_iterator histIter = h.begin(); histIter != h.end(); histIter++)
	{
		if(histIter->first > max)
		{
			max = histIter->first;
		}
		count += histIter->second;
		sum += (histIter->second * histIter->first);
	}
	
	double mean = sum / count;
	
	// Compute the standard devition
	sum = 0;
	count = 0;
	for(histogram::const_iterator histIter = h.begin(); histIter != h.end(); histIter++)
	{

		count += histIter->second;
		sum += (static_cast<double>(histIter->second) * pow((static_cast<double>(histIter->first) - mean), 2));
	}
	sum /= count;
	m_stdDev = sqrt(sum);	
	printf("mean distance: %lf std dev: %lf sum %lf count %lf\n", mean, m_stdDev, sum, count);
	

	// Create the initial pdf with all values being 0
	m_pdf = PDF(max+1, 0.0f);
	for(int i = 0; i < max+1; i++)
	{
		histogram::const_iterator iter = h.find(i);
		if(iter != h.end())
		{
			m_pdf[i] = static_cast<double>(iter->second) / count;
		}
		else
		{
			m_pdf[i] = MINP;
		}
	}
}

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
// Get the standard deviation of the estimate based on the empirical standard deviation and the number of data points
//
double Scaffold::GetStdDevOfEstimate(int n)
{
	return m_stdDev / sqrt((double)n);
}


//
//
//
range Scaffold::GenerateRange(int distance, int size, int n, int numDevs)
{
	double s = m_stdDev / double(n);
	int error = (int)s * numDevs;
	range r;
	r.start = distance - error;
	r.end = distance + error + size;
	return r;
}

PairAlign Scaffold::GetPairAlignsForRead(ReadID id)
{
	
	PairAlign pairAlign;
		
	ConstPMIter r2Iter = m_pairMap.find(id);
	ReadID id2;
	if(r2Iter != m_pairMap.end())
	{
		id2 = r2Iter->second;
	}
	else
	{
		// Not found (this shouldnt happen)
		assert(false);
	}
	
	// Get the read alignments
	ConstAMIter raIter = m_alignMap.find(id);
	assert(raIter != m_alignMap.end());
	pairAlign.pairs[0] = raIter->second;
	
	raIter = m_alignMap.find(id2);
	assert(raIter != m_alignMap.end());
	pairAlign.pairs[1] = raIter->second;
	
	return pairAlign;
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
	GenerateAllPairAlignments(contigID, cpvMap);
	
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
		GenerateAllPairAlignments(link.slaveID, cpvMap2);
		
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
		
		fileHandle >> id1 >> c1 >> pos1 >> rc1 >> id2 >> c2 >> pos2 >> rc2;

		// Check if this read was valid
		if(fileHandle.eof())
		{
			continue;
		}
		ReadAlign ra1 = BuildReadAlign(id1, c1, pos1, rc1);
		ReadAlign ra2 = BuildReadAlign(id2, c2, pos2, rc2);
		
		// Set up all the mappings
		m_alignMap[id1] = ra1;
		m_alignMap[id2] = ra2;
		
		m_pairMap[id1] = id2;
		m_pairMap[id2] = id1;
		
		m_contigReadMap[c1].push_back(id1);
		m_contigReadMap[c2].push_back(id2);
	}
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
		int size;
		Sequence seq;
		fileHandle >> head;
		
		if(head != '>' || fileHandle.eof())
		{
			continue;
		}		
		
		
		fileHandle >> contigID;
		fileHandle >> size;
		fileHandle >> seq;
		
		m_contigMap[contigID].seq = seq;
		m_contigMap[contigID].merged = false;
	}
	
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
// Print a linkage object
//
void Scaffold::PrintLinkage(ContigLinkage& link)
{
	std::string strOrient = (link.orientation == CORIEN_SAME) ? "same" : "opposite";
	std::string strOrder = (link.order == CORDER_LEFT) ? "left" : "right";
	cout << '\t' << link.slaveID << " orientation: " << strOrient << " c1 is on " << strOrder << " distance: " << link.distance << endl;	
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

