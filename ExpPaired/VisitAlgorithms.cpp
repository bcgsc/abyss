#include "VisitAlgorithms.h"
#include "Timer.h"
#include "SetOperations.h"

void ContigDataFunctions::merge(const ContigID& targetKey, ContigData& targetData, const ContigID& slaveKey, 
								ContigData& slaveData, extDirection dir, bool reverse, 
								bool removeChild, bool usableMerge)
{
	// Perform the merge
	targetData.merge(slaveData, dir, reverse, usableMerge);
	
	// Update the pairs on the master contig
	targetData.resolvePairs(m_pResolvePolicy, dir);
	targetData.resolvePairs(m_pResolvePolicy, !dir);
	
	// if the child is being removed, tell it to delete its alignments from the DB otherwise update its pairs
	if(removeChild)
	{
		slaveData.removeFromAlignDB(m_pDatabase);
	}
	else
	{
		slaveData.resolvePairs(m_pResolvePolicy, dir);
		slaveData.resolvePairs(m_pResolvePolicy, !dir);
	}
	
	(void)targetKey;
	(void)slaveKey;

	
	// validate the new target data
	targetData.validateSequences();	
	

}

void ContigDataFunctions::deleteCallback(const ContigID& slaveKey, const ContigData& slaveData)
{
	(void)slaveKey;
	(void)slaveData;
	printf("Implement the delete keys callback\n");
	assert(false);
}

bool ContigDataFunctions::check(ContigData& target, const ContigData& slave, extDirection dir, bool reverse)
{
	// should the slave be reversed?
	Sequence targetSeq = target.getSeq();
	Sequence slaveSeq = slave.getSeq();
	if(reverse)
	{
		slaveSeq = reverseComplement(slaveSeq);
	}
	
	Sequence* leftSeq;
	Sequence* rightSeq;
	// Order the contigs
	if(dir == SENSE)
	{
		leftSeq = &targetSeq;
		rightSeq = &slaveSeq;
	}
	else
	{
		leftSeq = &slaveSeq;
		rightSeq = &targetSeq;
	}
	
	// Get the last k bases of the left and the first k bases of the right
	PackedSeq leftEnd = leftSeq->substr(leftSeq->length() - m_overlap, m_overlap);
	PackedSeq rightBegin = rightSeq->substr(0, m_overlap);
	
	target.validateSequences();
	return leftEnd == rightBegin;
}

void DBGenerator::visit(const ContigID& /*id*/, const ContigData& data)
{
	// add a record for every sequence in the data
	data.writeToAlignDB(m_pDatabase);
}

void DBValidator::visit(const ContigID& /*id*/, const ContigData& data)
{
	// add a record for every sequence in the data
	data.validateDatabase(m_pDatabase);
}

void ContigDataOutputter::visit(const ContigID& id, const ContigData& data)
{
	// super hacky
	int iID = atoi(id.c_str());
	m_pWriter->WriteSequence(data.getSeq(), iID, 0);
}

int PairedMerger::resolve(DirectedGraph<ContigID, ContigData>* pGraph, const ContigID id)
{
	int numMerged = 0;
	printf("\n\n***ATTEMPTING TO RESOLVE %s***\n\n", id.c_str());
	
	// Functors
	SequenceDataCost dataCost(m_kmer);

	ContigDataFunctions dataFunctor(m_kmer, m_pResolvePolicy, m_pDatabase);
	
	for(size_t dirIdx = 1; dirIdx <= 1; ++dirIdx)
	{
		ContigIDSet reachableSet;
		extDirection dir = (extDirection)dirIdx;
		const ContigData& data = pGraph->getDataForVertex(id);
		
		// get the reachable set using the paired info 
		getReachableSet(data, dir, reachableSet);
		
		printf("Nodes reachable by pairs from %s: ", id.c_str());
		SetOperations::printSet(reachableSet);
		printf("\n");
		
		// Create the list of sequences that have been merged so far
		ContigIDSet mergedPathSet;
		
		// add the source node to the merged path set
		mergedPathSet.insert(id);
		
		// Get all the pairs described by this contig in the direction, even if they are from ambiguous contigs
		ContigSupportMap totalSupport;
		data.getSupportMap(dir, PSF_UNRESOLVED_ONLY, totalSupport);

		// enumerate all the possible paths that go through these nodes
		DirectedGraph<ContigID, ContigData>::FeasiblePaths possibleSolutions;
		
		// -1 is passed as the max number of paths as we would like everything
		pGraph->findSuperpaths(id, dir, reachableSet, possibleSolutions, dataCost, -1);
		
		// Score the paths
		DirectedGraph<ContigID, ContigData>::VertexPath chosenPath;
		bool found = false;
		size_t bestScore = 999;
		
		printf("Generated paths: \n");
		for(DirectedGraph<ContigID, ContigData>::FeasiblePaths::iterator solIter = possibleSolutions.begin(); solIter != possibleSolutions.end(); ++solIter)
		{
			size_t pathScore = pGraph->calculatePathLength(*solIter, dataCost);
			//size_t pathScore = scorePath(totalSupport, *solIter);
			printf("%zu ", pathScore);
			SetOperations::printPath(*solIter);
			printf("\n");
			
			if(pathScore < bestScore)
			{
				found = true;
				bestScore = pathScore;
				chosenPath = *solIter;
			}
		}
		assert(false && "clean this up to use the dirs passed back in the path");
		if(found)
		{
			printf("Chosen Path: \n");
			SetOperations::printPath(chosenPath);
			printf("\n");
			
			for(DirectedGraph<ContigID, ContigData>::VertexPath::const_iterator iter = chosenPath.begin(); iter != chosenPath.end(); ++iter)
			{
				// check if the node should be appended (target node stays) or merged (target node removed)
				ContigIDSet reachableSet;
				
				// Get the edge from the parent to the child
				EdgeDescription parentEdgeDesc = pGraph->getUniqueEdgeDesc(id, iter->key, dir);
				EdgeDescription relativeDesc = parentEdgeDesc.makeRelativeDesc();
								
				printf("Parent dir: %d relative dir: %d\n", parentEdgeDesc.dir, relativeDesc.dir);
				
				// Get the reachable set for the target node in the direction towards the parents
				ContigIDSet targetReachSet;
				ContigID currID = iter->key;
				const ContigData& currData = pGraph->getDataForVertex(currID);
				
				
				getReachableSet(currData, !relativeDesc.dir, targetReachSet);
				printf("Nodes reachable by pairs from %s: ", currID.c_str());
				SetOperations::printSet(targetReachSet);
				printf("\n");
				
				// Get the difference between the merged so far this set
				
				// If the resulting set is empty, the defined path is a superset of this node and the node can be removed
				ContigIDSet pathDiff;
				SetOperations::difference(targetReachSet, mergedPathSet, pathDiff);
				
				printf("Difference from merged path: ");
				SetOperations::printSet(pathDiff);
				printf("\n");
				
				// Count the number of times this node is in the path
				int pathMult = 0;
				for(DirectedGraph<ContigID, ContigData>::VertexPath::const_iterator countIter = chosenPath.begin(); countIter != chosenPath.end(); ++countIter)
				{
					if(countIter->key == currID)
					{
						pathMult++;
					}
				}
				
				// Decide what to do with the child contig
				// 1) If the pairs of this contig are completely resolved by the path leading to it should be removed
				// 2) If the pairs of this contig have not been resolved by any other path, it should be marked as usable (and considered to be equivalent in copy number to
				// the parent contig)
				
				// if the path child's reachable set is a strict subset of the parent's path, it can be removed
				bool removeChild = pathDiff.empty() && !targetReachSet.empty() && (pathMult == 1);
				
				// Check if the join between this sequence and the  is unique
				bool uniqueJoin = isJoinUnique(currData, !relativeDesc.dir, id);
				
				bool usableChild = removeChild && uniqueJoin;
				printf("%s remove? %d usable? %d\n", currID.c_str(), removeChild, usableChild);
				
				pGraph->mergePath(id, iter->key, dir, removeChild, usableChild, dataFunctor);
				
				mergedPathSet.insert(currID);
				numMerged++;
			}
		}
	}
	
	return numMerged;
}




//
// Generate the histogram of paired distances based on contigs that have self pairings
// We only use contigs greater than the distance cutoff or else the distribution will be biased towards smaller pairs
//
void PairedDistHistogramGenerator::visit(const ContigID& /*id*/, ContigData& data)
{
	(void)data;
	if(data.getLength() > m_distanceCutoff)
	{
		data.addSelfPairsToHist(m_pHist);
	}
}

//
// Set up the paired resolution policy, the criteria to decide upon whether a pair is resolved by a single contig or not
//
PairedResolvePolicy::PairedResolvePolicy(const PDF& pairedPdf)
{
	const double RESOLVED_PROP = 1.0f;
	// Calculate the minimal range that RESOLVED_PROP% of the pairs will fall into
	pairedPdf.calculateMinimalRange(RESOLVED_PROP, m_lowerLimit, m_upperLimit);
	
	printf("Resolution range is [%zu-%zu] for %lf of the pairs\n", m_lowerLimit, m_upperLimit, RESOLVED_PROP);
	
}



//
// Based on the policy, is this a resolved pairing
//
bool PairedResolvePolicy::isResolved(const ContigID& id1, ContigID& id2, size_t /*distance*/) const
{
	if(id1 == id2)// && distance >= m_lowerLimit && distance <= m_upperLimit)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//
// Functions not related to any particular functor
//




//
// Score a path based on the support map passed in
//
size_t scorePath(const ContigSupportMap& scoreMap, const ContigIDVec& path)
{
	size_t score = 0;
	for(ContigIDVec::const_iterator iter = path.begin(); iter != path.end(); ++iter)
	{
		// Look up the score for this key
		ContigSupportMap::const_iterator scoreIter = scoreMap.find(*iter);
		if(scoreIter != scoreMap.end())
		{
			score += scoreIter->second;
		}
	}
	
	score /= (path.size());
	return score;
}

//
// Get the set of reachable nodes (by pairs) from the specified contig in the specified direction
// Maybe move this to the contig data class
//
void getReachableSet(const ContigData& data, extDirection dir, ContigIDSet& idSet)
{
	// Get all the unresolved, usable pairs
	ContigSupportMap supportMap;
	data.getSupportMap(dir, PSF_USABLE_ONLY | PSF_UNRESOLVED_ONLY | PSF_UNIQUE_ALIGN_ONLY, supportMap);
	
	ContigData::supportMap2IDSet(supportMap, idSet);
	return;
}

//
// Check if the join between data and targetID satisfies all of contig data's pairs
//
bool isJoinUnique(const ContigData& data, extDirection dir, ContigID targetID)
{
	PairAlignmentsCollection pairAlignColl;
	data.extractAlignments(dir, PSF_ALL, pairAlignColl);

	ContigID selfID = data.getID();
	
	for(PairAlignmentsCollection::iterator colIter = pairAlignColl.begin(); colIter != pairAlignColl.end(); ++colIter)
	{
		bool passedAlign = false;
		// Ensure that every pair has an alignment pairing that is either (currID, currID) or (currID, targetID)
		for(PairAlignments::iterator apIter = colIter->begin(); apIter != colIter->end(); ++apIter)
		{
			ContigID pairContigID = apIter->pairSeqAlign.alignment.contigID;
			if(pairContigID == targetID || pairContigID == selfID)
			{
				passedAlign = true;
				break;
			}
		}
		
		if(!passedAlign)
		{
			return false;
		}
	}
	
	// If we reached here, all the contigs were resolved
	return true;
}

