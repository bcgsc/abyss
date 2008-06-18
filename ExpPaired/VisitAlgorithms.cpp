#include "VisitAlgorithms.h"
#include "Timer.h"
#include "SetOperations.h"

void ContigDataFunctions::merge(const ContigID& targetKey, ContigData& targetData, const ContigID& slaveKey, 
								ContigData& slaveData, extDirection dir, bool reverse, 
								bool usableMerge, bool removeChild)
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
	LinkFinder linkFinder;
	ContigDataFunctions dataFunctor(m_kmer, m_pResolvePolicy, m_pDatabase);
	
	for(size_t dirIdx = 0; dirIdx <= 1; ++dirIdx)
	{
		ContigIDSet reachableSet;
		extDirection dir = (extDirection)dirIdx;
		const ContigData& data = pGraph->getDataForVertex(id);
		
		// get the reachable set using the paired info 
		linkFinder.getReachableSet(id, data, dir, reachableSet);
		
		printf("Nodes reachable by pairs from %s: ", id.c_str());
		SetOperations::printSet(reachableSet);
		printf("\n");

		// check if (the shortest?) path from the root node to any node in the set is a super set of all the nodes
		ContigIDVec superPath;	
		bool found = pGraph->findSuperpath(id, dir, reachableSet, superPath, dataCost);

		// Create the list of sequences that have been merged so far
		ContigIDSet mergedPathSet;
		
		// add the source node to the merged path set
		mergedPathSet.insert(id);
		
		
		if(found)
		{
			for(ContigIDVec::const_iterator iter = superPath.begin(); iter != superPath.end(); ++iter)
			{
				// check if the node should be appended (target node stays) or merged (target node removed)
				ContigIDSet reachableSet;
				
				// Get the edge from the parent to the child
				EdgeDescription parentEdgeDesc = pGraph->getUniqueEdgeDesc(id, *iter, dir);
				EdgeDescription relativeDesc = parentEdgeDesc.makeRelativeDesc();
								
				printf("Parent dir: %d relative dir: %d\n", parentEdgeDesc.dir, relativeDesc.dir);
				
				// Get the reachable set for the target node in the direction towards the parents
				ContigIDSet targetReachSet;
				ContigID currID = *iter;
				const ContigData& currData = pGraph->getDataForVertex(currID);
				
				
				linkFinder.getReachableSet(currID, currData, !relativeDesc.dir, targetReachSet);
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
				for(ContigIDVec::iterator countIter = superPath.begin(); countIter != superPath.end(); ++countIter)
				{
					if(*countIter == currID)
					{
						pathMult++;
					}
				}
				
				// if the path child's reachable set is a strict subset of the parent's path, it can be removed
				bool removeChild = pathDiff.empty() && !targetReachSet.empty() && (pathMult == 1);
				//bool removeChild = (data.m_copyNumber == currData.m_copyNumber) && (pathMult == 1);

				pGraph->mergePath(id, *iter, dir, removeChild, dataFunctor);
				
				mergedPathSet.insert(currID);
				numMerged++;
			}
		}
	}
	
	return numMerged;
}


void LinkFinder::visit(const ContigID& /*id*/, ContigData& /*data*/)
{
	assert(false);
}

void LinkFinder::getReachableSet(const ContigID& /*id*/, const ContigData& data, extDirection dir, ContigIDSet& idSet)
{
	ContigSupportMap supportMap;
	data.getUnresolvedUsableSet(dir, supportMap);
	
	// Convert the support map to a set
	ContigIDSet reachableSet;
	for(ContigSupportMap::const_iterator iter = supportMap.begin(); iter != supportMap.end(); ++iter)
	{
		printf("Contig: %s Support: %d\n", iter->first.c_str(), iter->second);
		idSet.insert(iter->first);
	}
	
	return;
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

